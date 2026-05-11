"""
Configuration loading and validation for SEM calculations.
"""

import json
import sys
import logging
from pathlib import Path

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
except ImportError:
    rank = 0  # Assume serial execution if mpi4py not available

logger = logging.getLogger(__name__)

def load_config(config_path):
    """
    Load configuration from JSON file.
    
    Args:
        config_path: Path to JSON configuration file
        
    Returns:
        config: Parsed configuration dictionary
    """
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        if rank == 0:
            logger.info(f"Loaded configuration from {config_path}")
        return config
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
        if 'comm' in globals():
            comm.Abort(1)
        else:
            sys.exit(1)
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing JSON configuration: {e}")
        if 'comm' in globals():
            comm.Abort(1)
        else:
            sys.exit(1)

def validate_config(config, require_analyte=True):
    """
    Validate configuration dictionary.
    
    Args:
        config: Configuration dictionary
        require_analyte: If False, skip validation of moving analyte file existence
        
    Returns:
        bool: True if valid, False otherwise
    """
    required_sections = ["input", "pore_geometry", "simulation", "movement", "output"]
    
    for section in required_sections:
        if section not in config:
            logger.error(f"Missing required configuration section: {section}")
            return False
    
    # Validate pore type
    pore_type = config["pore_geometry"].get("pore_type", "").lower()
    valid_pore_types = ["cylindrical", "double_cone", "conical", "biological", "bin_file"]
    
    if pore_type not in valid_pore_types:
        logger.error(f"Invalid pore_type: {pore_type}. Must be one of: {valid_pore_types}")
        return False
    
    # Type-specific validation
    if pore_type == "biological" and not config["pore_geometry"].get("biological_pore_pdb"):
        logger.error("For biological pore, biological_pore_pdb must be provided")
        return False
    
    if pore_type == "bin_file" and not config["pore_geometry"].get("bin_file_path"):
        logger.error("For bin_file pore, bin_file_path must be provided")
        return False
    
    if pore_type == "double_cone":
        pore_radius = config["pore_geometry"].get("pore_radius", 100.0)
        outer_radius = config["pore_geometry"].get("outer_radius", pore_radius * 1.5)
        if outer_radius <= pore_radius:
            logger.error("For double_cone pore, outer_radius must be greater than pore_radius")
            return False

    if pore_type == "conical":
        top_radius = config["pore_geometry"].get("top_radius")
        bottom_radius = config["pore_geometry"].get("bottom_radius")
        if top_radius is None or bottom_radius is None:
            logger.error(
                "For conical pore, both top_radius and bottom_radius must be provided"
            )
            return False
        try:
            top_radius = float(top_radius)
            bottom_radius = float(bottom_radius)
        except (TypeError, ValueError):
            logger.error("conical top_radius and bottom_radius must be numeric")
            return False
        if top_radius <= 0 or bottom_radius <= 0:
            logger.error("conical top_radius and bottom_radius must be > 0")
            return False
    
    # Ensure defaults for optional parameters
    sim_section = config.get("simulation", {})
    if "xy_margin" not in sim_section:
        sim_section["xy_margin"] = 0.0
    else:
        try:
            sim_section["xy_margin"] = float(sim_section["xy_margin"])
        except (TypeError, ValueError):
            logger.error("Simulation parameter 'xy_margin' must be numeric")
            return False

    mesh_engine = str(sim_section.get("mesh_engine", "dolfinx")).lower()
    if mesh_engine not in ("dolfinx", "gmsh"):
        logger.error("Simulation parameter 'mesh_engine' must be 'dolfinx' or 'gmsh'")
        return False
    sim_section["mesh_engine"] = mesh_engine

    if mesh_engine == "gmsh":
        for key in ("gmsh_fine_size", "gmsh_coarse_size"):
            value = sim_section.get(key)
            if value is not None:
                try:
                    value = float(value)
                except (TypeError, ValueError):
                    logger.error("Simulation parameter '%s' must be numeric", key)
                    return False
                if value <= 0:
                    logger.error("Simulation parameter '%s' must be > 0", key)
                    return False
                sim_section[key] = value

        reproducible = sim_section.get("gmsh_reproducible", False)
        if isinstance(reproducible, (int, bool)):
            sim_section["gmsh_reproducible"] = bool(reproducible)
        else:
            logger.error("Simulation parameter 'gmsh_reproducible' must be a boolean")
            return False

        gmsh_num_threads = sim_section.get("gmsh_num_threads", None)
        if gmsh_num_threads is not None:
            try:
                gmsh_num_threads = int(gmsh_num_threads)
            except (TypeError, ValueError):
                logger.error("Simulation parameter 'gmsh_num_threads' must be an integer")
                return False
            if gmsh_num_threads <= 0:
                logger.error("Simulation parameter 'gmsh_num_threads' must be > 0")
                return False
            sim_section["gmsh_num_threads"] = gmsh_num_threads

        gmsh_random_seed = sim_section.get("gmsh_random_seed", None)
        if gmsh_random_seed is not None:
            try:
                gmsh_random_seed = int(gmsh_random_seed)
            except (TypeError, ValueError):
                logger.error("Simulation parameter 'gmsh_random_seed' must be an integer")
                return False
            if gmsh_random_seed < 0:
                logger.error("Simulation parameter 'gmsh_random_seed' must be >= 0")
                return False
            sim_section["gmsh_random_seed"] = gmsh_random_seed

        gmsh_random_factor = sim_section.get("gmsh_random_factor", None)
        if gmsh_random_factor is not None:
            try:
                gmsh_random_factor = float(gmsh_random_factor)
            except (TypeError, ValueError):
                logger.error("Simulation parameter 'gmsh_random_factor' must be numeric")
                return False
            if gmsh_random_factor < 0:
                logger.error("Simulation parameter 'gmsh_random_factor' must be >= 0")
                return False
            sim_section["gmsh_random_factor"] = gmsh_random_factor

        save_mesh_xdmf = sim_section.get("save_mesh_xdmf", False)
        if isinstance(save_mesh_xdmf, (int, bool)):
            sim_section["save_mesh_xdmf"] = bool(save_mesh_xdmf)
        else:
            logger.error("Simulation parameter 'save_mesh_xdmf' must be a boolean")
            return False

        fine_box = sim_section.get("gmsh_fine_box")
        if fine_box is not None:
            if not isinstance(fine_box, (list, tuple)) or len(fine_box) != 3:
                logger.error("Simulation parameter 'gmsh_fine_box' must be a 3-element list")
                return False
            try:
                sim_section["gmsh_fine_box"] = [float(v) for v in fine_box]
            except (TypeError, ValueError):
                logger.error("Simulation parameter 'gmsh_fine_box' must contain numeric values")
                return False

        fine_center = sim_section.get("gmsh_fine_center")
        if fine_center is not None:
            if not isinstance(fine_center, (list, tuple)) or len(fine_center) != 3:
                logger.error("Simulation parameter 'gmsh_fine_center' must be a 3-element list")
                return False
            try:
                sim_section["gmsh_fine_center"] = [float(v) for v in fine_center]
            except (TypeError, ValueError):
                logger.error("Simulation parameter 'gmsh_fine_center' must contain numeric values")
                return False

        fine_center_mode = str(sim_section.get("gmsh_fine_center_mode", "origin")).lower()
        if fine_center_mode not in ("origin", "analyte_com"):
            logger.error("Simulation parameter 'gmsh_fine_center_mode' must be 'origin' or 'analyte_com'")
            return False
        sim_section["gmsh_fine_center_mode"] = fine_center_mode

    # Validate input files exist
    input_pdb = config["input"]["moving_pdb"]
    if require_analyte:
        if not Path(input_pdb).exists():
            logger.error(f"Input PDB file not found: {input_pdb}")
            return False
    
    if pore_type == "biological":
        bio_pdb = config["pore_geometry"]["biological_pore_pdb"]
        if bio_pdb and not Path(bio_pdb).exists():
            logger.error(f"Biological pore PDB file not found: {bio_pdb}")
            return False
    
    if pore_type == "bin_file":
        bin_file = config["pore_geometry"]["bin_file_path"]
        if bin_file and not Path(bin_file).exists():
            logger.error(f"Binary file not found: {bin_file}")
            return False
    
    if rank == 0:
        logger.info("Configuration validation passed")
    return True

def print_config_summary(config):
    """
    Print a summary of the loaded configuration.
    Enhanced to show all pore types including binary files.
    
    Args:
        config: Configuration dictionary
    """
    if rank == 0:
        logger.info("Configuration Summary:")
        logger.info(f"  Input PDB: {config['input']['moving_pdb']}")
        
        pore_geom = config["pore_geometry"]
        logger.info(f"  Pore Type: {pore_geom['pore_type']}")
        
        if pore_geom['pore_type'].lower() == "biological":
            logger.info(f"  Biological Pore PDB: {pore_geom.get('biological_pore_pdb', 'Not specified')}")
        elif pore_geom['pore_type'].lower() == "bin_file":
            logger.info(f"  Binary File: {pore_geom.get('bin_file_path', 'Not specified')}")
            logger.info(f"  Binary file stores: {pore_geom.get('bin_file_units', 'distance')}")
            if "mask_radius" in pore_geom and pore_geom["mask_radius"] > 0:
                logger.info(f"  Mask Radius: {pore_geom['mask_radius']} Å")
        elif pore_geom['pore_type'].lower() == "conical":
            logger.info(f"  Top Radius: {pore_geom.get('top_radius', 'Not specified')} Å")
            logger.info(f"  Bottom Radius: {pore_geom.get('bottom_radius', 'Not specified')} Å")
        else:
            logger.info(f"  Pore Radius: {pore_geom.get('pore_radius', 100.0)} Å")
            if pore_geom.get("corner_radius", 0) > 0:
                logger.info(f"  Corner Radius: {pore_geom['corner_radius']} Å")
            if "outer_radius" in pore_geom:
                logger.info(f"  Outer Radius: {pore_geom['outer_radius']} Å")
        
        logger.info(f"  Membrane Thickness: {pore_geom['membrane_thickness']} Å")
        
        sim = config["simulation"]
        logger.info(f"  Voltage: {sim['voltage']} mV")
        logger.info(f"  Bulk Conductivity: {sim['bulk_conductivity']} S/m")
        logger.info(f"  Use VdW Radii: {sim['use_vdw_radii']}")
        logger.info(f"  Use internal PQR prep: {sim['use_pdb2pqr']}")
        logger.info(f"  XY Margin: {sim.get('xy_margin', 0.0)} Å")
        mesh_engine = sim.get("mesh_engine", "dolfinx")
        logger.info(f"  Mesh Engine: {mesh_engine}")
        if str(mesh_engine).lower() == "gmsh":
            logger.info(f"  Gmsh Fine Size: {sim.get('gmsh_fine_size', 'default')} Å")
            logger.info(f"  Gmsh Coarse Size: {sim.get('gmsh_coarse_size', 'default')} Å")
            logger.info(f"  Gmsh Fine Box: {sim.get('gmsh_fine_box', None)}")
            logger.info(f"  Gmsh Fine Center: {sim.get('gmsh_fine_center', None)}")
            logger.info(f"  Gmsh Fine Center Mode: {sim.get('gmsh_fine_center_mode', 'origin')}")
            logger.info(f"  Gmsh Reproducible: {sim.get('gmsh_reproducible', False)}")
            logger.info(f"  Gmsh Threads: {sim.get('gmsh_num_threads', 'default')}")
            logger.info(f"  Gmsh Random Seed: {sim.get('gmsh_random_seed', 'default')}")
            logger.info(f"  Gmsh Random Factor: {sim.get('gmsh_random_factor', 'default')}")
        logger.info(f"  Save Mesh XDMF: {sim.get('save_mesh_xdmf', False)}")
        if "membrane_conductivity" in sim:
            logger.info(f"  Membrane Conductivity: {sim['membrane_conductivity']} S/m")
        if sim.get("prevent_analyte_overlap"):
            logger.info("  Analyte overlap protection: enabled")
            if sim.get("use_radius_overlap_check", False):
                logger.info(
                    "  Radius-sum overlap check: enabled (buffer %.3f Å)",
                    float(sim.get("overlap_buffer", 0.0)),
                )
            if sim.get("overlap_distance_threshold") is not None:
                logger.info(
                    "  Fixed-distance overlap threshold: %.3f Å",
                    float(sim.get("overlap_distance_threshold")),
                )
        if "cleanup_temp_files" in sim:
            logger.info(f"  Cleanup temporary files: {sim['cleanup_temp_files']}")

        movement = config["movement"]
        logger.info(f"  Z Range: {movement['z_start']} to {movement['z_end']} Å")
        logger.info(f"  Z Step: {movement['z_step']} Å")
        
        output = config["output"]
        logger.info(f"  Output Prefix: {output['output_prefix']}")
        logger.info(f"  Preview Frames: {output.get('preview_frames', 0)}")
        if "preview_mesh" in output:
            logger.info(f"  Preview Mesh: {output['preview_mesh']}")
        if "arbd_export" in output and output["arbd_export"]:
            arbd = output["arbd_export"]
            ions = arbd.get("ions", [])
            logger.info(
                "  ARBD export: ions=%s, stride=%s, wall=%s kcal/mol, T=%s K",
                ions,
                arbd.get("stride", 0),
                arbd.get("wall_height", 100.0),
                arbd.get("temperature_K", 295.0),
            )

def create_example_config(pore_type="cylindrical", output_file="example_config.json"):
    """
    Create an example configuration file.
    
    Args:
        pore_type: Type of pore to create example for
        output_file: Output filename
    """
    base_config = {
        "input": {
            "moving_pdb": "analyte.pdb"
        },
        "simulation": {
            "voltage": 500.0,
            "bulk_conductivity": 1.660843,
            "grid_resolution": 1.0,
            "use_pdb2pqr": True,  # Use internal PDBFixer-based PQR preparation
            "use_vdw_radii": False,
            "default_radius": 1.5,
            "membrane_conductivity": 0.0001,
            "cleanup_temp_files": True,
            "use_charges": False,  # Toggle charge-aware modulation
            "debye_length": 2.15,  # Angstroms, for ~2M KCl
            "bjerrum_length": 7.15,  # Angstroms, for water at 293K
            "charge_clip": 2.0,  # Clip cosh factor to [1/clip, clip] e.g. [0.5, 2.0]
            "xy_margin": 0.0,  # Additional XY padding applied to mesh extent (Å)
            "mesh_engine": "dolfinx",
            "gmsh_reproducible": False,
            "gmsh_num_threads": None,
            "gmsh_random_seed": None,
            "gmsh_random_factor": None,
            "save_mesh_xdmf": False,
            "prevent_analyte_overlap": False,
            "use_radius_overlap_check": False,
            "overlap_buffer": 0.0,
            "overlap_distance_threshold": None,
        },
        "movement": {
            "z_start": 150.0,
            "z_end": -150.0,
            "z_step": 1.0
        },
        "output": {
            "output_prefix": f"{pore_type}_sem",
            "preview_frames": 4
        }
    }
    
    if pore_type == "cylindrical":
        base_config["pore_geometry"] = {
            "pore_type": "cylindrical",
            "pore_radius": 100.0,
            "corner_radius": 0.0,
            "membrane_thickness": 200.0
        }
    elif pore_type == "double_cone":
        base_config["pore_geometry"] = {
            "pore_type": "double_cone",
            "pore_radius": 80.0,
            "outer_radius": 120.0,
            "membrane_thickness": 200.0
        }
    elif pore_type == "conical":
        base_config["pore_geometry"] = {
            "pore_type": "conical",
            "top_radius": 120.0,
            "bottom_radius": 60.0,
            "membrane_thickness": 200.0
        }
    elif pore_type == "biological":
        base_config["pore_geometry"] = {
            "pore_type": "biological",
            "biological_pore_pdb": "pore_structure.pdb",
            "membrane_thickness": 200.0
        }
    elif pore_type == "bin_file":
        base_config["pore_geometry"] = {
            "pore_type": "bin_file",
            "bin_file_path": "pore_structure.bin",
            "mask_radius": -1,
            "membrane_thickness": 200.0,
            "bin_file_units": "distance"
        }
    
    with open(output_file, 'w') as f:
        json.dump(base_config, f, indent=2)
    
    if rank == 0:
        logger.info(f"Created example {pore_type} configuration: {output_file}")
