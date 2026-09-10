"""
Command line interface for SEM calculations.
"""

import argparse
import copy
import csv
import json
import logging
import sys
from pathlib import Path
from typing import Optional

from . import __version__ as SEM_VERSION
from .config import load_config, validate_config, print_config_summary, create_example_config
from .grid_io import RAMP_MAX, RAMP_MIN
from .provenance import git_commit as _git_commit
# `VerticalMovementSEM`/`AnalyteOverlapError` transitively import dolfinx
# (mpi4py, petsc4py, ufl, ...), which only ships via conda-forge. Importing
# them at module scope would break `import sem.cli` (and therefore
# `create_sem_from_config`'s argument-parsing-only callers, e.g.
# `create_config`) in a pip-only environment that never installs dolfinx.
# Imported lazily where actually needed: inside `create_sem_from_config`
# and `run_rotation_scan` below.
from .rotation import (
    rotate_pdb_to_grid_center,
    parse_angle_file,
    random_uniform_rotations,
    rotation_matrix_from_spec,
    RotationSpec,
)

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
except ImportError:
    comm = None
    rank = 0  # Assume serial execution if mpi4py not available

logger = logging.getLogger(__name__)


def _open_pore_result_lines(sem, open_current):
    """
    Build the full `{prefix}_open_pore_current.txt` file content as a list
    of lines (every line but the last is a `#`-prefixed comment; the last
    line is the bare current value). Pure string formatting over `sem`'s
    public attributes -- no dolfinx import -- so it is unit-testable with a
    lightweight stand-in object instead of a real (dolfinx-backed)
    `VerticalMovementSEM` instance. Kept in sync with the console summary
    printed alongside it in `main()`.
    """
    lines = ["# Open pore current calculation results", f"# Pore type: {sem.pore_type}"]

    if sem.pore_type == "cylindrical":
        lines.append(f"# Pore radius: {sem.pore_radius:.1f} Angstrom")
        if sem.corner_radius > 0:
            lines.append(f"# Corner radius: {sem.corner_radius:.1f} Angstrom")
            chamfer_depth = sem.chamfer_depth if sem.chamfer_depth is not None else sem.corner_radius
            lines.append(f"# Chamfer depth: {chamfer_depth:.1f} Angstrom")
    elif sem.pore_type == "double_cone":
        lines.append(f"# Inner radius: {sem.pore_radius:.1f} Angstrom")
        lines.append(f"# Outer radius: {sem.outer_radius:.1f} Angstrom")
    elif sem.pore_type == "conical":
        lines.append(f"# Top radius: {sem.top_radius:.1f} Angstrom")
        lines.append(f"# Bottom radius: {sem.bottom_radius:.1f} Angstrom")
    elif sem.pore_type == "profile":
        lines.append(f"# Profile: {sem.profile_path}")
    elif sem.pore_type == "biological":
        lines.append(f"# Biological pore: {sem.biological_pore_pdb}")
    elif sem.pore_type == "bin_file":
        lines.append(f"# Binary file: {sem.bin_file_path}")
        lines.append(f"# Binary file units: {sem.bin_file_units}")

    if sem.pore_type in ("cylindrical", "double_cone", "conical", "profile"):
        lines.append(f"# Distance metric: {sem.distance_metric}")

    lines.append(f"# Membrane thickness: {sem.membrane_thickness:.1f} Angstrom")
    lines.append(f"# Applied voltage: {sem.voltage*1000:.1f} mV")
    lines.append(f"# Bulk conductivity: {sem.bulk_conductivity:.1f} S/m")
    lines.append(f"# Grid resolution: {sem.grid_resolution:.1f} Angstrom")
    lines.append(f"# sem_version: {SEM_VERSION}")
    lines.append(f"# git_commit: {_git_commit()}")
    lines.append(f"# ramp: {RAMP_MIN} {RAMP_MAX}")
    lines.append("# Open_pore_current(nA)")
    lines.append(f"{open_current:.6e}")
    return lines


def create_sem_from_config(config, prepare_analyte=True, *, gmsh_center_mode_override=None):
    """
    Create VerticalMovementSEM instance from configuration dictionary.
    Enhanced to support all pore types including binary files.
    
    Args:
        config: Configuration dictionary loaded from JSON
        prepare_analyte: Whether to load and prepare the moving analyte structure
        
    Returns:
        sem: VerticalMovementSEM instance
        config: Configuration dictionary (returned for convenience)
    """
    # `VerticalMovementSEM` transitively imports dolfinx; deferred here so
    # `sem.cli` stays importable (and `create_config` usable) without it.
    from .vertical_movement_sem import VerticalMovementSEM

    if rank == 0:
        logger.info("Creating SEM instance from configuration...")

    # Extract parameters from config
    moving_pdb = config["input"]["moving_pdb"]

    # Pore geometry parameters
    pore_geom = config["pore_geometry"]
    pore_type = pore_geom["pore_type"].lower()  # Convert to lowercase

    # Common parameters
    pore_radius = pore_geom.get("pore_radius", 100.0)
    membrane_thickness = pore_geom["membrane_thickness"]

    # Type-specific parameters
    corner_radius = pore_geom.get("corner_radius", 0.0)
    chamfer_depth = pore_geom.get("chamfer_depth", None)
    distance_metric = pore_geom.get("distance_metric", "euclidean")
    profile_path = pore_geom.get("profile_path", None)
    outer_radius = pore_geom.get("outer_radius", None)
    top_radius = pore_geom.get("top_radius", None)
    bottom_radius = pore_geom.get("bottom_radius", None)
    biological_pore_pdb = pore_geom.get("biological_pore_pdb", None)
    bin_file_path = pore_geom.get("bin_file_path", None)
    bin_file_units = pore_geom.get("bin_file_units", "distance").lower()
    mask_radius = pore_geom.get("mask_radius", -1)
    membrane_z_offset = pore_geom.get("membrane_z_offset")
    
    # Simulation parameters
    sim = config["simulation"]
    voltage = sim["voltage"]
    bulk_conductivity = sim["bulk_conductivity"]
    grid_resolution = sim["grid_resolution"]
    use_vdw_radii = sim["use_vdw_radii"]
    use_pdb2pqr = sim.get("use_pdb2pqr")
    use_charges = sim.get("use_charges", False)
    debye_length = sim.get("debye_length", 2.15)
    bjerrum_length = sim.get("bjerrum_length", 7.15)
    charge_clip = sim.get("charge_clip", 2.0)
    default_radius = sim["default_radius"]
    xy_margin = float(sim.get("xy_margin", 0.0))
    prevent_analyte_overlap = bool(sim.get("prevent_analyte_overlap", False))
    use_radius_overlap_check = bool(sim.get("use_radius_overlap_check", False))
    overlap_buffer = float(sim.get("overlap_buffer", 0.0))
    overlap_distance_threshold = sim.get("overlap_distance_threshold", None)
    mesh_engine = sim.get("mesh_engine", "dolfinx")
    gmsh_fine_size = sim.get("gmsh_fine_size", None)
    gmsh_coarse_size = sim.get("gmsh_coarse_size", None)
    gmsh_fine_box = sim.get("gmsh_fine_box", None)
    gmsh_fine_center = sim.get("gmsh_fine_center", None)
    gmsh_fine_center_mode = sim.get("gmsh_fine_center_mode", "origin")
    gmsh_reproducible = sim.get("gmsh_reproducible", False)
    gmsh_num_threads = sim.get("gmsh_num_threads", None)
    gmsh_random_seed = sim.get("gmsh_random_seed", None)
    gmsh_random_factor = sim.get("gmsh_random_factor", None)
    save_mesh_xdmf = sim.get("save_mesh_xdmf", False)
    if gmsh_center_mode_override is not None:
        gmsh_fine_center_mode = gmsh_center_mode_override
        if gmsh_center_mode_override == "origin":
            gmsh_fine_center = None
    
    # Optional parameters
    membrane_conductivity = sim.get("membrane_conductivity", 0.0001)
    cleanup_temp_files = sim.get("cleanup_temp_files", True)
    

    # Movement parameters
    movement = config["movement"]
    z_start = movement["z_start"]
    z_end = movement["z_end"]
    z_step = movement["z_step"]
    
    # Output parameters
    output = config["output"]
    output_prefix = output["output_prefix"]
    # Optional ARBD-compatible DX export block; passed through to VerticalMovementSEM.
    arbd_export_cfg = output.get("arbd_export", None)
    if arbd_export_cfg and rank == 0:
        ions = arbd_export_cfg.get("ions", [])
        logger.info(
            "ARBD export enabled: ions=%s, stride=%s, wall_height=%s kcal/mol, T=%s K",
            ions,
            arbd_export_cfg.get("stride", 0),
            arbd_export_cfg.get("wall_height", 100.0),
            arbd_export_cfg.get("temperature_K", 295.0),
        )
    
    # Box dimensions (optional)
    box_dims = config.get("box_dimensions", None)
    if box_dims:
        box_dimensions = {
            'x': tuple(box_dims["x"]),
            'y': tuple(box_dims["y"]),
            'z': tuple(box_dims["z"])
        }
        if rank == 0:
            logger.info("Using box dimensions from configuration")
    else:
        box_dimensions = None
        if rank == 0:
            logger.info("Box dimensions will be auto-calculated")
    
    # Create SEM instance
    sem = VerticalMovementSEM(
        moving_pdb=moving_pdb,
        pore_type=pore_type,
        pore_radius=pore_radius,
        outer_radius=outer_radius,
        top_radius=top_radius,
        bottom_radius=bottom_radius,
        corner_radius=corner_radius,
        chamfer_depth=chamfer_depth,
        distance_metric=distance_metric,
        profile_path=profile_path,
        biological_pore_pdb=biological_pore_pdb,
        bin_file_path=bin_file_path,
        mask_radius=mask_radius,
        membrane_thickness=membrane_thickness,
        z_start=z_start,
        z_end=z_end,
        z_step=z_step,
        voltage=voltage,
        output_prefix=output_prefix,
        box_dimensions=box_dimensions,
        grid_resolution=grid_resolution,
        bulk_conductivity=bulk_conductivity,
        use_vdw_radii=use_vdw_radii,
        use_pdb2pqr=use_pdb2pqr,
        use_charges=use_charges,
        debye_length=debye_length,
        bjerrum_length=bjerrum_length,
        charge_clip=charge_clip,
        default_radius=default_radius,
        membrane_conductivity=membrane_conductivity,
        membrane_z_offset=membrane_z_offset,
        cleanup_temp_files=cleanup_temp_files,
        xy_margin=xy_margin,
        mesh_engine=mesh_engine,
        gmsh_fine_size=gmsh_fine_size,
        gmsh_coarse_size=gmsh_coarse_size,
        gmsh_fine_box=gmsh_fine_box,
        gmsh_fine_center=gmsh_fine_center,
        gmsh_fine_center_mode=gmsh_fine_center_mode,
        gmsh_reproducible=gmsh_reproducible,
        gmsh_num_threads=gmsh_num_threads,
        gmsh_random_seed=gmsh_random_seed,
        gmsh_random_factor=gmsh_random_factor,
        save_mesh_xdmf=save_mesh_xdmf,
        prepare_analyte=prepare_analyte,
        prevent_analyte_overlap=prevent_analyte_overlap,
        use_radius_overlap_check=use_radius_overlap_check,
        overlap_buffer=overlap_buffer,
        overlap_distance_threshold=overlap_distance_threshold,
        bin_file_units=bin_file_units,
        arbd_export=arbd_export_cfg,
    )
    
    if rank == 0:
        logger.info("SEM instance created successfully")
    return sem, config


def _resolve_path(path_str: Optional[str], base_dir: Path) -> Optional[Path]:
    if path_str is None:
        return None
    candidate = Path(path_str)
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()
    return candidate


def _broadcast(obj):
    try:
        return comm.bcast(obj, root=0)  # type: ignore[name-defined]
    except NameError:
        return obj


def run_rotation_scan(base_config: dict, args: argparse.Namespace, config_file: Path | None):
    # Deferred: see the comment on the top-level import removal above.
    from .vertical_movement_sem import AnalyteOverlapError

    if rank == 0:
        logger.info("Starting rotation scan with mode '%s'", args.mode)

    base_dir = config_file.parent if config_file else Path.cwd()
    base_config_copy = copy.deepcopy(base_config)

    base_moving_pdb = _resolve_path(base_config_copy["input"]["moving_pdb"], base_dir)
    dx_path = _resolve_path(getattr(args, "dx", None), base_dir)
    angles_path = _resolve_path(getattr(args, "angles", None), base_dir)

    if base_moving_pdb is None:
        raise ValueError("Configuration is missing input.moving_pdb")

    if rank == 0:
        if not base_moving_pdb.exists():
            raise FileNotFoundError(f"Input PDB not found: {base_moving_pdb}")
        if dx_path is not None and not dx_path.exists():
            raise FileNotFoundError(f"DX file not found: {dx_path}")
        if angles_path and not angles_path.exists():
            raise FileNotFoundError(f"Angles file not found: {angles_path}")

    if comm is not None and rank != 0:
        rotations = None
    else:
        if angles_path:
            rotations = list(parse_angle_file(angles_path))
            if args.start_index:
                rotations = rotations[args.start_index:]
            if args.samples is not None:
                rotations = rotations[: args.samples]
        else:
            if args.samples is None:
                raise ValueError("rotation_scan requires --samples when no angles file is provided.")
            rotations = random_uniform_rotations(args.samples, seed=args.seed)
            if args.start_index:
                rotations = rotations[args.start_index:]

    if comm is not None:
        rotations = _broadcast(rotations)

    reuse_mesh = getattr(args, "reuse_mesh", False)
    shared_sem_instance = None
    reuse_open_pore = getattr(args, "reuse_open_pore", False)
    shared_open_current = None

    base_prepare_analyte = args.mode != 'open_pore'
    if reuse_mesh:
        if rank == 0 and args.mode == 'open_pore':
            logger.warning("reuse-mesh is enabled, but open_pore mode ignores analyte orientation.")
        base_prepare_analyte = True
        gmsh_center_override = None
        if args.mode in ("open_pore", "preview_only"):
            gmsh_center_override = "origin"
        shared_sem_instance, _ = create_sem_from_config(
            base_config_copy,
            prepare_analyte=base_prepare_analyte,
            gmsh_center_mode_override=gmsh_center_override,
        )

    if args.mode == "run" and reuse_open_pore:
        if rank == 0:
            logger.info("Computing open pore current once for rotation scan.")
        if reuse_mesh and shared_sem_instance is not None:
            shared_open_current = shared_sem_instance.calculate_open_pore_current()
        else:
            sem_for_open, _ = create_sem_from_config(
                base_config_copy,
                prepare_analyte=base_prepare_analyte,
            )
            shared_open_current = sem_for_open.calculate_open_pore_current()
        if rank == 0:
            logger.info("Open pore current cached: %.6e nA", shared_open_current)

    if not rotations:
        if rank == 0:
            logger.warning("No rotations available (check angles file or sample limits).")
        return

    output_dir = Path(args.output_dir).resolve()
    if rank == 0:
        output_dir.mkdir(parents=True, exist_ok=True)

    results: list[tuple[int, float, float, float, float]] = []
    # Hybrid trace: rows of (idx, rx, ry, rz, z, current, normalized_current, blockage),
    # populated only when mode == 'run' so PETK can render a (z × rotation) overlay.
    hybrid_rows: list[tuple[int, float, float, float, float, float, float, float]] = []
    prefix_base = base_config_copy["output"].get("output_prefix", "vertical_movement")

    for offset, rotation_spec in enumerate(rotations):
        idx = args.start_index + offset
        rotation_dir = output_dir / f"rot_{idx:03d}"
        output_prefix_path = rotation_dir / f"{prefix_base}_rot_{idx:03d}"

        if rank == 0:
            rotation_dir.mkdir(parents=True, exist_ok=True)

        if reuse_mesh:
            if rank == 0:
                config_variant = copy.deepcopy(base_config_copy)
                config_variant.setdefault("rotation", {})
                config_variant["rotation"].update(
                    {"rx": rotation_spec.rx, "ry": rotation_spec.ry, "rz": rotation_spec.rz}
                )
                config_variant["output"]["output_prefix"] = str(output_prefix_path)
                config_path = rotation_dir / "config.json"
                with open(config_path, "w") as handle:
                    json.dump(config_variant, handle, indent=2)
            else:
                config_variant = None
        else:
            if rank == 0:
                rotated_name = f"moving_{rotation_spec.label()}.pdb"
                rotated_path = rotation_dir / rotated_name
                rotate_pdb_to_grid_center(
                    base_moving_pdb,
                    dx_path,
                    rotation_spec,
                    rotated_path,
                )

                config_variant = copy.deepcopy(base_config_copy)
                config_variant["input"]["moving_pdb"] = str(rotated_path)
                config_variant["output"]["output_prefix"] = str(output_prefix_path)
                config_variant.setdefault("rotation", {})
                config_variant["rotation"].update(
                    {"rx": rotation_spec.rx, "ry": rotation_spec.ry, "rz": rotation_spec.rz}
                )

                config_path = rotation_dir / "config.json"
                with open(config_path, "w") as handle:
                    json.dump(config_variant, handle, indent=2)
            else:
                config_variant = None

        if comm is not None:
            config_variant = _broadcast(config_variant)
            comm.Barrier()

        if reuse_mesh:
            sem_instance = shared_sem_instance
            sem_instance.output_prefix = str(output_prefix_path)
            rotation_matrix = rotation_matrix_from_spec(rotation_spec)
            sem_instance.set_analyte_rotation_matrix(rotation_matrix)
        else:
            gmsh_center_override = None
            if args.mode in ("open_pore", "preview_only"):
                gmsh_center_override = "origin"
            sem_instance, _ = create_sem_from_config(
                config_variant,
                prepare_analyte=base_prepare_analyte,
                gmsh_center_mode_override=gmsh_center_override,
            )

        try:
            if args.mode == 'open_pore':
                open_current = sem_instance.calculate_open_pore_current()
                if rank == 0:
                    results.append((idx, rotation_spec.rx, rotation_spec.ry, rotation_spec.rz, open_current))
            elif args.mode == 'run':
                if shared_open_current is not None:
                    run_results = sem_instance.run(open_current=shared_open_current)
                else:
                    run_results = sem_instance.run()
                if rank == 0 and run_results and len(run_results.get("currents", [])) > 0:
                    final_current = float(run_results["currents"][-1])
                    results.append((idx, rotation_spec.rx, rotation_spec.ry, rotation_spec.rz, final_current))
                    # Capture full z-trace for hybrid output. Lengths are aligned by run().
                    z_positions = run_results.get("z_positions")
                    if z_positions is None:
                        z_positions = []
                    currents_full = run_results.get("currents")
                    if currents_full is None:
                        currents_full = []
                    norm_full = run_results.get("normalized_currents")
                    if norm_full is None:
                        norm_full = []
                    block_full = run_results.get("blockages")
                    if block_full is None:
                        block_full = []
                    n_z = len(z_positions)
                    if n_z and n_z == len(currents_full):
                        # Pad missing optional series with NaN so columns align.
                        if len(norm_full) != n_z:
                            norm_full = [float("nan")] * n_z
                        if len(block_full) != n_z:
                            block_full = [float("nan")] * n_z
                        for k in range(n_z):
                            hybrid_rows.append((
                                idx,
                                float(rotation_spec.rx),
                                float(rotation_spec.ry),
                                float(rotation_spec.rz),
                                float(z_positions[k]),
                                float(currents_full[k]),
                                float(norm_full[k]),
                                float(block_full[k]),
                            ))
            elif args.mode == 'preview_only':
                if rank == 0:
                    logger.info("Generating preview frames for rotation %s", rotation_spec.label())
                preview_frames = config_variant["output"].get("preview_frames", 4)
                from .visualization import create_preview_frames, export_dx_file, export_mesh
                create_preview_frames(sem_instance, num_frames=preview_frames)
                if "preview_dx" in config_variant["output"]:
                    export_dx_file(sem_instance, config_variant["output"]["preview_dx"])
                if "preview_mesh" in config_variant["output"]:
                    mesh_path = config_variant["output"]["preview_mesh"]
                    if isinstance(mesh_path, bool):
                        mesh_path = None
                    export_mesh(sem_instance, mesh_path)
            else:
                raise ValueError(f"Unsupported rotation mode: {args.mode}")
        except AnalyteOverlapError as overlap_exc:
            if rank == 0:
                logger.warning(
                    "Skipping rotation %s due to analyte overlap: %s",
                    rotation_spec.label(),
                    overlap_exc,
                )
            continue

    if reuse_mesh and shared_sem_instance is not None:
        shared_sem_instance.reset_analyte_rotation()

    if results and rank == 0:
        results_path = output_dir / "rotation_results.csv"
        with open(results_path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["index", "rx", "ry", "rz", "current_nA"])
            writer.writerows(results)
        logger.info("Rotation scan results saved to %s", results_path)

    if hybrid_rows and rank == 0:
        hybrid_path = output_dir / "hybrid_currents.csv"
        with open(hybrid_path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([
                "index", "rx", "ry", "rz",
                "z_A", "current_nA", "normalized_current", "blockage",
            ])
            writer.writerows(hybrid_rows)
        logger.info(
            "Hybrid (translocation × rotation) trace saved to %s "
            "(%d rows across %d rotations)",
            hybrid_path, len(hybrid_rows), len(results),
        )

def main(argv=None):
    """
    Main function to run SEM with JSON configuration.
    Supports multiple modes: 'run', 'preview_only', 'open_pore', and 'create_config'
    Enhanced to support all pore types including binary files.

    Args:
        argv: Optional list of command-line arguments (excluding the program
            name), e.g. ``["create_config", "conical", "-o", "cfg.json"]``.
            Defaults to ``sys.argv[1:]`` when omitted (normal CLI usage).
    """
    parser = argparse.ArgumentParser(
        description='Run SEM calculation with JSON configuration (supports all pore types including binary files)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m sem config.json run           # Run full simulation
  python -m sem config.json preview_only  # Generate preview plots only
  python -m sem config.json open_pore     # Calculate open pore current only
  python -m sem config.json rotation_scan map.dx angles.txt --samples 10
  python -m sem create_config cylindrical # Create example config file
  python -m sem derive_geometry ...       # Derive geometry from an all-atom map

Pore Types:
  - cylindrical: Simple cylindrical pore with optional corner rounding
  - double_cone: Hourglass-shaped pore
  - conical: Single frustum (asymmetric cone) pore
  - profile: Pore wall defined by an arbitrary radius-vs-z table (CSV)
  - biological: PDB-based biological pore structure
  - bin_file: Binary file-based pore structure (like original code)
        """
    )
    
    # Add subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Run/preview/open_pore commands
    run_parser = subparsers.add_parser('run', help='Run SEM simulation')
    run_parser.add_argument('config', help='Path to JSON configuration file')
    
    preview_parser = subparsers.add_parser('preview_only', help='Generate preview plots only')
    preview_parser.add_argument('config', help='Path to JSON configuration file')
    
    open_pore_parser = subparsers.add_parser('open_pore', help='Calculate open pore current only')
    open_pore_parser.add_argument('config', help='Path to JSON configuration file')

    rotation_parser = subparsers.add_parser('rotation_scan', help='Rotate analyte and run SEM for each orientation')
    rotation_parser.add_argument('config', help='Base JSON configuration file')
    rotation_parser.add_argument('dx', nargs='?', help='DX file describing the reference grid used for rotation pivot (optional)')
    rotation_parser.add_argument('angles', nargs='?',
                                 help='Text file containing rx ry rz (degrees) per line; if omitted, random sampling is used')
    rotation_parser.add_argument('--samples', type=int, default=None,
                                 help='Maximum number of rotations to process (or required when angles file omitted)')
    rotation_parser.add_argument('--output-dir', default='rotations',
                                 help='Directory to store rotated structures and results (default: rotations)')
    rotation_parser.add_argument('--mode', choices=['open_pore', 'run', 'preview_only'], default='open_pore',
                                 help='SEM mode to execute for each rotation (default: open_pore)')
    rotation_parser.add_argument('--start-index', type=int, default=0,
                                 help='Index offset used in naming rotations (default: 0)')
    rotation_parser.add_argument('--seed', type=int, default=None,
                                 help='Random seed used when generating orientations without an angles file')
    rotation_parser.add_argument('--reuse-mesh', action='store_true',
                                 help='Rotate analyte in-place and reuse a single SEM mesh instance (no PDB rewriting)')
    rotation_parser.add_argument('--reuse-open-pore', action='store_true',
                                 help='Compute open pore current once and reuse for all rotations (assumes mesh is unchanged)')
    
    # Create config command
    config_parser = subparsers.add_parser('create_config', help='Create example configuration file')
    config_parser.add_argument(
        'pore_type',
        choices=['cylindrical', 'double_cone', 'conical', 'biological', 'bin_file', 'profile'],
        help='Type of pore for example configuration')
    config_parser.add_argument('-o', '--output', default='example_config.json',
                              help='Output filename (default: example_config.json)')

    # Derive-geometry command: fits parametric pore parameters (or a
    # radius-vs-z profile table) from an all-atom distance/conductivity map.
    # The implementation lives in sem.scripts.derive_geometry, imported
    # lazily (see below) so `sem.cli` itself stays importable without
    # scipy; this subparser exists purely for discoverability (--help).
    # Actual dispatch bypasses argparse below so its own argument parsing
    # (unknown to this module) is untouched.
    subparsers.add_parser(
        'derive_geometry',
        help='Derive parametric pore parameters (or a profile table) from an all-atom map',
    )

    raw_argv = list(sys.argv[1:]) if argv is None else list(argv)

    # `derive_geometry` forwards its remaining arguments verbatim to
    # sem.scripts.derive_geometry.main, whose own CLI is defined in that
    # module -- so it is dispatched here, before argparse gets a chance to
    # interpret those arguments itself. The import is lazy (done here
    # rather than at module scope) so that `sem.cli` stays importable in
    # an environment without scipy, which only `derive_geometry` needs.
    if raw_argv and raw_argv[0] == 'derive_geometry':
        try:
            from .scripts.derive_geometry import main as derive_geometry_main
        except ImportError as exc:
            print(
                f"Error: 'derive_geometry' is not available: {exc}",
                file=sys.stderr,
            )
            sys.exit(2)
        return derive_geometry_main(raw_argv[1:])

    # Legacy support for old command line format
    if len(raw_argv) >= 2 and raw_argv[1] in ['run', 'preview_only', 'open_pore']:
        # Old format: program config.json run
        args = parser.parse_args([raw_argv[1], raw_argv[0]])
    else:
        args = parser.parse_args(raw_argv)

    if args.command == 'create_config':
        create_example_config(args.pore_type, args.output)
        return

    if args.command is None:
        if rank == 0:
            parser.print_help()
        return
    
    config_path: Path | None = Path(args.config).resolve() if hasattr(args, "config") else None
    
    try:
        # Load and validate configuration
        config = load_config(args.config)
        
        require_analyte = args.command != 'open_pore'

        if not validate_config(config, require_analyte=require_analyte):
            if rank == 0:
                logger.error("Configuration validation failed")
            sys.exit(1)
        
        # Print configuration summary
        print_config_summary(config)

        if args.command == 'rotation_scan':
            run_rotation_scan(config, args, config_path)
            if rank == 0:
                logger.info("Execution completed successfully!")
            return
        
        # Create SEM instance
        prepare_analyte = args.command != 'open_pore'
        gmsh_center_override = None
        if args.command in ("open_pore", "preview_only"):
            gmsh_center_override = "origin"
        sem, config = create_sem_from_config(
            config,
            prepare_analyte=prepare_analyte,
            gmsh_center_mode_override=gmsh_center_override,
        )
        
        # Run based on mode
        if args.command == 'open_pore':
            if rank == 0:
                logger.info("Running in open pore current calculation mode")
            
            # Calculate open pore current only
            try:
                open_current = sem.calculate_open_pore_current()
                
                # Print results to console
                if rank == 0:
                    print(f"\n{'='*60}")
                    print(f"OPEN PORE CURRENT CALCULATION RESULTS")
                    print(f"{'='*60}")
                    print(f"Pore type: {sem.pore_type}")
                    if sem.pore_type == "cylindrical":
                        print(f"Pore radius: {sem.pore_radius:.1f} Å")
                        if sem.corner_radius > 0:
                            print(f"Corner radius: {sem.corner_radius:.1f} Å")
                            chamfer_depth = sem.chamfer_depth if sem.chamfer_depth is not None else sem.corner_radius
                            print(f"Chamfer depth: {chamfer_depth:.1f} Å")
                    elif sem.pore_type == "double_cone":
                        print(f"Inner radius: {sem.pore_radius:.1f} Å")
                        print(f"Outer radius: {sem.outer_radius:.1f} Å")
                    elif sem.pore_type == "conical":
                        print(f"Top radius: {sem.top_radius:.1f} Å")
                        print(f"Bottom radius: {sem.bottom_radius:.1f} Å")
                    elif sem.pore_type == "profile":
                        print(f"Profile: {sem.profile_path}")
                    elif sem.pore_type == "biological":
                        print(f"Biological pore: {sem.biological_pore_pdb}")
                    elif sem.pore_type == "bin_file":
                        print(f"Binary file: {sem.bin_file_path}")
                        print(f"Binary file units: {sem.bin_file_units}")

                    if sem.pore_type in ("cylindrical", "double_cone", "conical", "profile"):
                        print(f"Distance metric: {sem.distance_metric}")

                    print(f"Membrane thickness: {sem.membrane_thickness:.1f} Å")
                    print(f"Applied voltage: {sem.voltage*1000:.1f} mV")
                    print(f"Bulk conductivity: {sem.bulk_conductivity:.1f} S/m")
                    print(f"Grid resolution: {sem.grid_resolution:.1f} Å")
                    print(f"Open pore current: {open_current:.6e} nA")
                    print(f"SEM version: {SEM_VERSION}")
                    print(f"Git commit: {_git_commit()}")
                    print(f"Ramp: {RAMP_MIN}-{RAMP_MAX} Å")
                    print(f"{'='*60}")

                    # Save results to file
                    output_file = f"{sem.output_prefix}_open_pore_current.txt"
                    with open(output_file, 'w') as f:
                        f.write("\n".join(_open_pore_result_lines(sem, open_current)) + "\n")

                    logger.info(f"Results saved to: {output_file}")
                    logger.info("Open pore current calculation completed successfully!")
                
            except Exception as e:
                if rank == 0:
                    logger.error(f"Error calculating open pore current: {e}")
                    import traceback
                    logger.error(f"Traceback: {traceback.format_exc()}")
                sys.exit(1)
        
        elif args.command == 'preview_only':
            if rank == 0:
                logger.info("Running in preview-only mode")
            preview_frames = config["output"].get("preview_frames", 4)
            
            # Import visualization functions
            try:
                from .visualization import create_preview_frames, export_dx_file, export_mesh
                
                # Create preview frames (now works for all pore types)
                create_preview_frames(sem, num_frames=preview_frames, save_plots=True)
                
                # Try to export DX file (now works for all pore types including bin_file)
                try:
                    dx_file = export_dx_file(sem)
                    if dx_file and rank == 0:
                        logger.info(f"DX file created: {dx_file}")
                        logger.info("Load in VMD with: vmd " + dx_file)
                except Exception as e:
                    if rank == 0:
                        logger.warning(f"Could not create DX file: {e}")

                preview_mesh = config["output"].get("preview_mesh", None)
                if preview_mesh is not None:
                    if preview_mesh is True:
                        preview_mesh = None
                    try:
                        mesh_file = export_mesh(sem, preview_mesh)
                        if mesh_file and rank == 0:
                            logger.info("Mesh file created: %s", mesh_file)
                    except Exception as e:
                        if rank == 0:
                            logger.warning("Could not create mesh file: %s", e)
                
                if rank == 0:
                    logger.info("Preview mode completed successfully!")
                
            except ImportError as e:
                if rank == 0:
                    logger.error(f"Could not import visualization functions: {e}")
                    logger.error("Make sure matplotlib is installed for preview functionality")
            
        elif args.command == 'run':
            if rank == 0:
                logger.info("Running full simulation")
            # Run full simulation
            results = sem.run()
        
        if rank == 0:
            logger.info("Execution completed successfully!")
        
    except Exception as e:
        if rank == 0:
            logger.error(f"Error during execution: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
        sys.exit(1)

if __name__ == "__main__":
    main()
