"""
Main SEM calculation class for vertical movement of analyte through nanopore.
Enhanced version with van der Waals radii support and multiple pore geometries.
Converted to DOLFINx.
"""

import numpy as np
import logging
import time
from typing import Optional
from typing import Optional
import dolfinx
import dolfinx.fem as fem
import dolfinx.mesh as mesh
import dolfinx.io as io
from dolfinx.fem.petsc import LinearProblem
import ufl
from mpi4py import MPI
from petsc4py import PETSc
import MDAnalysis as mda
from scipy.spatial import KDTree
from .structure_preparation import prepare_structure, PreparedStructure

from .utils import loadFunc, get_dof_coordinates
from .van_der_waals import VanDerWaalsRadii
from .pore_geometry import PoreGeometry
from .conductivity_models import SimpleConductivityModel, ChargeAwareConductivityModel

logger = logging.getLogger(__name__)


class AnalyteOverlapError(RuntimeError):
    """Raised when an analyte atom overlaps with membrane/pore walls."""


def _log_element_radius_summary(
    atoms,
    radii,
    *,
    logger_obj,
    prefix: str = "Analyte",
    detailed: bool = True,
):
    """
    Log a compact summary of radii used for each element in the provided atom group.
    """
    if logger_obj is None:
        return

    summary = {}
    for atom, radius in zip(atoms, radii):
        element = ""
        if hasattr(atom, "element") and atom.element:
            element = atom.element.strip()
        if not element and hasattr(atom, "name") and atom.name:
            element = atom.name.strip()
        element_key = element.upper() if element else "UNKNOWN"

        entry = summary.setdefault(element_key, {"count": 0, "radii": set()})
        entry["count"] += 1
        entry["radii"].add(round(float(radius), 4))

    logger_obj.info("%s radii summary:", prefix)

    if detailed:
        for element_key in sorted(summary):
            entry = summary[element_key]
            radii_list = ", ".join(f"{value:.4f}" for value in sorted(entry["radii"]))
            logger_obj.info(
                "  %s: %d atoms, radii %s",
                element_key,
                entry["count"],
                radii_list,
            )
    else:
        # Aggregate by primary element symbol (first letter, or original key if unknown)
        aggregate: dict[str, dict[str, set[float]]] = {}
        for element_key, entry in summary.items():
            primary = element_key[0] if element_key and element_key[0].isalpha() else element_key
            primary = primary.upper()
            agg_entry = aggregate.setdefault(primary, {"count": 0, "radii": set()})
            agg_entry["count"] += entry["count"]
            agg_entry["radii"].update(entry["radii"])

        for primary in sorted(aggregate):
            entry = aggregate[primary]
            radii_list = ", ".join(f"{value:.4f}" for value in sorted(entry["radii"]))
            logger_obj.info(
                "  %s: %d atoms, radii %s",
                primary,
                entry["count"],
                radii_list,
            )

class VerticalMovementSEM:
    """
    SEM calculation for vertical movement of analyte through nanopore.
    Enhanced version with van der Waals radii support and multiple pore geometries including biological pores and binary files.
    All units in Angstroms.
    Converted to DOLFINx.
    
    This class:
    1. Creates a mathematical conductivity grid for membrane with pore (cylindrical, double cone, biological, or binary file)
    2. Loads a moving analyte structure
    3. Assigns van der Waals radii based on atomic elements
    4. Generates virtual trajectory by moving analyte vertically
    5. Calculates conductance at each position
    6. Solves FEM for current at each position
    """
    
    def __init__(self, 
                 moving_pdb,        # Centered analyte
                 pore_type="cylindrical",  # "cylindrical", "double_cone", "conical", "biological", or "bin_file"
                 pore_radius=100.0,  # Pore radius (Å) - for cylindrical or inner radius for double cone
                 outer_radius=None,  # Outer radius for double cone (Å) - if None, uses pore_radius * 1.5
                 top_radius=None,    # Top-face radius for conical pore (Å)
                 bottom_radius=None, # Bottom-face radius for conical pore (Å)
                 corner_radius=0.0,  # Corner radius for cylindrical pore (Å)
                 biological_pore_pdb=None,  # Path to PDB file for biological pore
                 bin_file_path=None,  # Path to binary file for bin_file pore
                 bin_file_units="distance",  # How to interpret bin file values
                 mask_radius=-1,  # Mask radius for bin_file pore
                 membrane_thickness=200.0,  # Membrane thickness (Å)
                 z_start=150.0,     # Starting Z position (Å)
                 z_end=-150.0,      # Ending Z position (Å)
                 z_step=1.0,        # Step size (Å)
                 voltage=100.0,     # Applied voltage (mV)
                 output_prefix="vertical_movement",
                 box_dimensions=None,  # If None, auto-calculate
                 grid_resolution=1.0,  # Grid resolution (Å)
                 cutoff=5.0,          # Distance cutoff for analyte (Å)
                 bulk_conductivity=11.2,  # Bulk conductivity (S/m)
                 use_vdw_radii=True,  # Use van der Waals radii
                 use_pdb2pqr=False,
                 use_charges=False,  # Toggle charge-aware modulation
                 debye_length=2.15,  # Angstroms, for ~2M KCl
                 bjerrum_length=7.15,  # Angstroms, for water at 293K
                 charge_clip=2.0,  # Clip cosh factor to [1/clip, clip] e.g. [0.5, 2.0]
                 force_field="CHARMM",
                 ph=7.0,
                 default_radius=1.5,  # Default radius for unknown elements (Å)
                 membrane_conductivity=1e-7*11.2,
                 membrane_z_offset=0.0,
                 xy_margin=0.0,
                 mesh_engine="dolfinx",
                 gmsh_fine_size=None,
                 gmsh_coarse_size=None,
                 gmsh_fine_box=None,
                 gmsh_fine_center=None,
                 gmsh_fine_center_mode="origin",
                 gmsh_reproducible=False,
                 gmsh_num_threads=None,
                 gmsh_random_seed=None,
                 gmsh_random_factor=None,
                 save_mesh_xdmf=False,
                 cleanup_temp_files=True,
                 prepare_analyte=True,
                 prevent_analyte_overlap=False,
                 use_radius_overlap_check=False,
                 overlap_buffer=0.0,
                 overlap_distance_threshold=None):  # Overlap buffer (Å)
        
        # Initialize MPI
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        
        self.moving_pdb = moving_pdb
        self.pore_type = pore_type.lower()
        self.pore_radius = pore_radius
        self.outer_radius = outer_radius if outer_radius is not None else pore_radius * 1.5
        self.top_radius = top_radius
        self.bottom_radius = bottom_radius
        self.corner_radius = corner_radius
        self.biological_pore_pdb = biological_pore_pdb
        self.bin_file_path = bin_file_path
        self.bin_file_units = bin_file_units.lower() if bin_file_units else "distance"
        self.mask_radius = mask_radius
        self.membrane_thickness = membrane_thickness
        self.z_start = z_start
        self.z_end = z_end
        self.z_step = z_step
        self.voltage = voltage / 1000.0  # Convert to V
        self.output_prefix = output_prefix
        self.grid_resolution = grid_resolution
        self.cutoff = cutoff
        self.input_bulk_conductivity = float(bulk_conductivity)
        self.bulk_scaling_factor = 0.1
        self.bulk_conductivity = self.input_bulk_conductivity * self.bulk_scaling_factor
        self.use_vdw_radii = use_vdw_radii
        self.use_pdb2pqr = use_pdb2pqr
        self.force_field = force_field
        self.ph = ph
        self.use_charges = use_charges
        self.debye_length = debye_length
        self.bjerrum_length = bjerrum_length
        self.charge_clip = charge_clip
        self.base_phi_interp = None  # Will store pore potential interpolator if needed
        self.base_dist_interp = None  # Optional distance interpolator (bin_file)
        self.default_radius = default_radius
        self.membrane_conductivity = self.bulk_scaling_factor * membrane_conductivity
        self.membrane_z_offset = membrane_z_offset
        self.xy_margin = float(xy_margin)
        self.mesh_engine = (mesh_engine or "dolfinx").lower()
        if self.mesh_engine not in ("dolfinx", "gmsh"):
            raise ValueError("mesh_engine must be 'dolfinx' or 'gmsh'")
        self.gmsh_fine_size = gmsh_fine_size
        self.gmsh_coarse_size = gmsh_coarse_size
        self.gmsh_fine_box = None
        if gmsh_fine_box is not None:
            if len(gmsh_fine_box) != 3:
                raise ValueError("gmsh_fine_box must be a 3-element list [hx, hy, hz]")
            self.gmsh_fine_box = [float(v) for v in gmsh_fine_box]
        self.gmsh_fine_center = None
        if gmsh_fine_center is not None:
            if len(gmsh_fine_center) != 3:
                raise ValueError("gmsh_fine_center must be a 3-element list [x, y, z]")
            self.gmsh_fine_center = [float(v) for v in gmsh_fine_center]
        self.gmsh_fine_center_mode = (gmsh_fine_center_mode or "origin").lower()
        if self.gmsh_fine_center_mode not in ("origin", "analyte_com"):
            raise ValueError("gmsh_fine_center_mode must be 'origin' or 'analyte_com'")
        self.gmsh_reproducible = bool(gmsh_reproducible)
        self.gmsh_num_threads = int(gmsh_num_threads) if gmsh_num_threads is not None else None
        self.gmsh_random_seed = int(gmsh_random_seed) if gmsh_random_seed is not None else None
        self.gmsh_random_factor = float(gmsh_random_factor) if gmsh_random_factor is not None else None
        if self.gmsh_reproducible:
            if self.gmsh_num_threads is None:
                self.gmsh_num_threads = 1
            if self.gmsh_random_seed is None:
                self.gmsh_random_seed = 1
            if self.gmsh_random_factor is None:
                self.gmsh_random_factor = 0.0
        self.save_mesh_xdmf = bool(save_mesh_xdmf)
        self._gmsh_fine_center_user = self.gmsh_fine_center
        self.update_mesh_each_step = (
            self.mesh_engine == "gmsh"
            and self.gmsh_fine_center_mode == "analyte_com"
            and self._gmsh_fine_center_user is None
        )
        self._gmsh_mesh_center = None
        self._gmsh_center_warned = False
        self._gmsh_option_warned = set()
        self._mesh_write_count = 0
        self.cleanup_temp_files = cleanup_temp_files
        self.verbose_output = not self.cleanup_temp_files
        self.prepare_analyte = prepare_analyte
        self.analyte_prepared = False
        self.prevent_analyte_overlap = bool(prevent_analyte_overlap)
        self.use_radius_overlap_check = bool(use_radius_overlap_check)
        self.overlap_buffer = float(overlap_buffer)
        self.overlap_distance_threshold = (
            float(overlap_distance_threshold)
            if overlap_distance_threshold is not None
            else None
        )

        # Initialize bin file attributes
        self.bin_dimensions = None
        self.bin_grid_shape = None
        self.bin_grid_spacing = None
        if self.rank == 0:
            logger.info(
                "Bulk conductivity scaled by fixed factor %.6f (input %.4f S/m → %.4f S/m)",
                self.bulk_scaling_factor,
                self.input_bulk_conductivity,
                self.bulk_conductivity,
            )
            if self.mesh_engine == "gmsh":
                fine = self.gmsh_fine_size if self.gmsh_fine_size is not None else self.grid_resolution
                coarse = self.gmsh_coarse_size if self.gmsh_coarse_size is not None else max(self.grid_resolution * 3.0, self.grid_resolution)
                logger.info(
                    "Using Gmsh mesh (fine %.3f Å, coarse %.3f Å, fine box %s, fine center %s)",
                    fine,
                    coarse,
                    self.gmsh_fine_box,
                    self.gmsh_fine_center if self.gmsh_fine_center is not None else self.gmsh_fine_center_mode,
                )
                if self.gmsh_reproducible:
                    logger.info(
                        "Gmsh reproducibility enabled (threads=%s, seed=%s, random_factor=%s).",
                        self.gmsh_num_threads,
                        self.gmsh_random_seed,
                        self.gmsh_random_factor,
                    )
                if self.save_mesh_xdmf:
                    logger.info("Mesh XDMF output enabled.")
                if self.update_mesh_each_step:
                    logger.info("Gmsh fine center follows analyte COM; rebuilding mesh each position.")
            if self.prevent_analyte_overlap:
                logger.info("Analyte overlap protection enabled.")
                if self.use_radius_overlap_check:
                    logger.info(
                        "Radius-sum overlap check enabled (buffer %.3f Å).",
                        self.overlap_buffer,
                    )
                if self.overlap_distance_threshold is not None:
                    logger.info(
                        "Fixed-distance overlap threshold: %.3f Å.",
                        self.overlap_distance_threshold,
                    )
            if not self.verbose_output:
                logger.info("Per-position logging disabled (cleanup_temp_files=True). "
                            "Set cleanup_temp_files=False for detailed output dumps.")
        
        # Validate pore type / bin units
        if self.pore_type not in ["cylindrical", "double_cone", "conical", "biological", "bin_file"]:
            raise ValueError("pore_type must be 'cylindrical', 'double_cone', 'conical', 'biological', or 'bin_file'")
        if self.pore_type == "bin_file":
            valid_units = ("distance", "conductivity")
            if self.bin_file_units not in valid_units:
                raise ValueError(f"bin_file_units must be one of {valid_units}, got {self.bin_file_units!r}")
            if self.rank == 0:
                logger.info("Interpreting binary grid values as %s", self.bin_file_units)
        
        # For biological pore, ensure PDB file is provided
        if self.pore_type == "biological" and not self.biological_pore_pdb:
            raise ValueError("For biological pore, biological_pore_pdb must be provided")
        
        # For bin_file pore, ensure bin file is provided
        if self.pore_type == "bin_file" and not self.bin_file_path:
            raise ValueError("For bin_file pore, bin_file_path must be provided")
        
        # For double cone, ensure outer_radius > pore_radius (inner_radius)
        if self.pore_type == "double_cone" and self.outer_radius <= self.pore_radius:
            raise ValueError("For double_cone pore, outer_radius must be greater than pore_radius (inner_radius)")
        
        self.moving_universe = None
        self.moving_positions = np.zeros((0, 3))
        self.moving_radii = np.zeros(0)
        self.moving_charges = np.zeros(0)
        self.moving_com = np.zeros(3)
        self._base_moving_positions = None
        self._current_rotation_matrix = np.eye(3)
        self._open_pore_current = None

        if self.prepare_analyte:
            if self.rank == 0:
                logger.info("Loading analyte structure...")

            prepared_analyte: Optional[PreparedStructure] = None
            moving_file = moving_pdb
            try:
                if self.use_pdb2pqr:
                    if self.rank == 0:
                        logger.info(
                            "Preparing analyte with external pdb2pqr pipeline "
                            "and custom radius overrides."
                        )
                    try:
                        prepared_analyte = prepare_structure(
                            moving_pdb,
                            ph=self.ph,
                            default_radius=self.default_radius,
                            use_external_pdb2pqr=True,
                            pdb2pqr_force_field=self.force_field,
                        )
                    except ImportError as exc:
                        if self.rank == 0:
                            logger.error(
                                "Structure preparation requires pdbfixer and openmm: %s", exc
                            )
                        raise
                    moving_file = prepared_analyte.pqr_file

                self.moving_universe = mda.Universe(str(moving_file))
                moving_atoms = self.moving_universe.atoms
                self.moving_positions = moving_atoms.positions.copy()
                self.moving_com = moving_atoms.center_of_mass()
                self._base_moving_positions = self.moving_positions.copy()
                if self.rank == 0:
                    logger.info(f"Loaded analyte with {len(self.moving_positions)} atoms")

                if self.use_pdb2pqr:
                    try:
                        self.moving_radii = moving_atoms.radii.copy()
                        if self.rank == 0:
                            logger.info("Using radii from prepared PQR file")
                    except Exception as exc:
                        if self.rank == 0:
                            logger.error(f"Failed to read radii from prepared PQR: {exc}")
                        raise
                elif self.use_vdw_radii:
                    if self.rank == 0:
                        logger.info("Assigning van der Waals radii to analyte atoms...")
                    self.moving_radii = VanDerWaalsRadii.assign_radii_to_atoms(
                        moving_atoms,
                        default_radius=self.default_radius,
                        verbose=(self.rank == 0)
                    )
                else:
                    try:
                        self.moving_radii = moving_atoms.radii
                        if self.rank == 0:
                            logger.info("Using radii from analyte PDB file")
                    except Exception:
                        self.moving_radii = np.ones(len(moving_atoms)) * self.default_radius
                        if self.rank == 0:
                            logger.info(f"Using default radius {self.default_radius} Å for all analyte atoms")

                if self.rank == 0:
                    try:
                        _log_element_radius_summary(
                            moving_atoms,
                            self.moving_radii,
                            logger_obj=logger,
                            prefix="Analyte",
                            detailed=not self.cleanup_temp_files,
                        )
                    except Exception as exc:
                        logger.debug("Failed to log analyte radii summary: %s", exc)

                if self.use_pdb2pqr:
                    self.moving_charges = moving_atoms.charges.copy()
                    if self.use_charges and self.rank == 0 and not np.any(self.moving_charges):
                        logger.warning(
                            "Prepared PQR contains zero charges; charge-aware conductivity "
                            "model will use neutral atoms."
                        )
                else:
                    try:
                        self.moving_charges = moving_atoms.charges.copy()
                    except Exception:
                        self.moving_charges = np.zeros(len(moving_atoms))

                if self.rank == 0:
                    logger.info(f"Loaded analyte with {np.sum(self.moving_charges != 0)} charged atoms")
                self.analyte_prepared = True
            except Exception as exc:
                if self.rank == 0:
                    logger.error(f"Failed to load analyte file {moving_file}: {exc}")
                raise
            finally:
                if prepared_analyte is not None:
                    prepared_analyte.cleanup()
        else:
            if self.rank == 0:
                logger.info("Skipping analyte preparation (open pore mode)")
        # Set or calculate box dimensions
        if box_dimensions is None:
            self.calculate_box_dimensions()
        else:
            self.box_dimensions = box_dimensions
            
        # Create base conductivity grid for membrane
        self.create_base_conductivity_grid()
        
        # Initialize conductivity model for analyte
        if self.use_charges and self.prepare_analyte:
            self.conductivity_model = ChargeAwareConductivityModel(
                bulk_conductivity=self.bulk_conductivity,
                cutoff=self.cutoff,
                charge_clip=self.charge_clip
            )
        else:
            self.conductivity_model = SimpleConductivityModel(
                bulk_conductivity=self.bulk_conductivity,
            )
        
        # Setup DOLFINx
        self.setup_dolfinx()
        
    def calculate_box_dimensions(self):
        """Auto-calculate box dimensions based on pore geometry."""
        # Use fixed dimensions that encompass the membrane and movement range
        padding = 20.0  # Å
        
        # XY dimensions based on pore size and padding
        if self.pore_type == "cylindrical":
            max_radius = self.pore_radius
        elif self.pore_type == "double_cone":
            max_radius = self.outer_radius
        elif self.pore_type == "conical":
            if self.top_radius is None or self.bottom_radius is None:
                raise ValueError(
                    "Conical pore requires both top_radius and bottom_radius "
                    "to auto-calculate box dimensions."
                )
            max_radius = max(self.top_radius, self.bottom_radius)
        elif self.pore_type == "bin_file":
            # For bin files, try to read the dimensions
            try:
                from .utils import readbinGrid
                val3d, [Lm, Wm, Hm], [nx, ny, nz], metadata = readbinGrid(
                    self.bin_file_path, return_metadata=True
                )
                # Store for later use
                self.bin_dimensions = [Lm, Wm, Hm]
                self.bin_grid_shape = [nx, ny, nz]
                if metadata and "spacing" in metadata:
                    self.bin_grid_spacing = float(metadata["spacing"][0])
                max_radius = max(Lm, Wm) / 2 + padding
                if self.rank == 0:
                    logger.info(f"Bin file dimensions: {Lm:.1f} x {Wm:.1f} x {Hm:.1f} Å")
            except:
                max_radius = 100.0  # Default if reading fails
                if self.rank == 0:
                    logger.warning("Could not read bin file dimensions, using default")
        else:  # biological
            # For biological pores, use a reasonable default and adjust based on pore structure if needed
            max_radius = 50.0  # Default radius, can be adjusted
            if self.biological_pore_pdb:
                try:
                    # Quick check of pore dimensions
                    temp_universe = mda.Universe(self.biological_pore_pdb)
                    temp_positions = temp_universe.atoms.positions
                    temp_com = temp_universe.atoms.center_of_mass()
                    temp_positions -= temp_com  # Center
                    
                    # Calculate approximate radius from centered positions
                    max_dist = np.max(np.sqrt(temp_positions[:, 0]**2 + temp_positions[:, 1]**2))
                    max_radius = max(max_radius, max_dist + 20.0)  # Add some padding
                    if self.rank == 0:
                        logger.info(f"Estimated biological pore radius: {max_dist:.1f} Å")
                except:
                    if self.rank == 0:
                        logger.warning("Could not estimate biological pore size, using default")
            
        xy_size = max(150.0, max_radius * 3 + padding)
        
        # Z dimension based on movement range and membrane
        z_min = min(-self.membrane_thickness/2 - padding, self.z_end - padding)
        z_max = max(self.membrane_thickness/2 + padding, self.z_start + padding)
        
        self.box_dimensions = {
            'x': (-xy_size, xy_size),
            'y': (-xy_size, xy_size),
            'z': (z_min, z_max)
        }
        
        if self.rank == 0:
            logger.info(f"Box dimensions: X={self.box_dimensions['x']}, "
                       f"Y={self.box_dimensions['y']}, Z={self.box_dimensions['z']}")
    
    def create_base_conductivity_grid(self):
        if self.rank == 0:
            logger.info(f"Creating base conductivity grid for {self.pore_type} pore...")
        
        # Create grid if needed (for grid-based pores)
        if self.pore_type in ["cylindrical", "double_cone", "conical", "biological"]:
            x_range = np.linspace(
                self.box_dimensions['x'][0],
                self.box_dimensions['x'][1],
                int(round((self.box_dimensions['x'][1] - self.box_dimensions['x'][0]) / self.grid_resolution)) + 1
            )
            y_range = np.linspace(
                self.box_dimensions['y'][0],
                self.box_dimensions['y'][1],
                int(round((self.box_dimensions['y'][1] - self.box_dimensions['y'][0]) / self.grid_resolution)) + 1
            )
            z_range = np.linspace(
                self.box_dimensions['z'][0],
                self.box_dimensions['z'][1],
                int(round((self.box_dimensions['z'][1] - self.box_dimensions['z'][0]) / self.grid_resolution)) + 1
            )
            X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing='ij')
        else:
            X = Y = Z = None
        
        # Create pore object
        pore_kwargs = {}
        if self.pore_type in ["cylindrical", "double_cone", "conical", "biological"]:
            pore_kwargs['bulk_conductivity'] = self.bulk_conductivity
            pore_kwargs['membrane_conductivity'] = self.membrane_conductivity

        if self.pore_type == "cylindrical":
            pore_kwargs.update({
                'pore_radius': self.pore_radius,
                'membrane_half_thickness': self.membrane_thickness / 2,
                'corner_radius': self.corner_radius,
            })
        elif self.pore_type == "double_cone":
            pore_kwargs.update({
                'inner_radius': self.pore_radius,
                'outer_radius': self.outer_radius,
                'membrane_half_thickness': self.membrane_thickness / 2
            })
        elif self.pore_type == "conical":
            if self.top_radius is None or self.bottom_radius is None:
                raise ValueError(
                    "Conical pore requires both top_radius and bottom_radius."
                )
            pore_kwargs.update({
                'top_radius': self.top_radius,
                'bottom_radius': self.bottom_radius,
                'membrane_half_thickness': self.membrane_thickness / 2
            })
        elif self.pore_type == "bin_file":
            pore_kwargs.update({
                'bin_file_path': self.bin_file_path,
                'base_sigma': self.bulk_conductivity,
                'mask_radius': self.mask_radius,
                'data_units': self.bin_file_units,
            })
        elif self.pore_type == "biological":
            pore_kwargs.update({
                'pore_pdb': self.biological_pore_pdb,
                'membrane_half_thickness': self.membrane_thickness / 2,
                'cutoff': self.cutoff,
                'use_vdw_radii': self.use_vdw_radii,
                'default_radius': self.default_radius,
                'membrane_z_offset': self.membrane_z_offset,
                'use_charges': self.use_charges,
                'debye_length': self.debye_length,
                'bjerrum_length': self.bjerrum_length,
                'resolution': self.grid_resolution,
                'cleanup_temp_files': self.cleanup_temp_files,
                'box_dimensions': self.box_dimensions,
                'use_direct_distance_calculation': False,  # Or your default
                'use_pdb2pqr': self.use_pdb2pqr,
                'force_field': self.force_field,
                'ph': self.ph
            })
        
        self.pore_obj = PoreGeometry.create_pore(self.pore_type, X=X, Y=Y, Z=Z, **pore_kwargs)
        self.base_cond_interp = self.pore_obj.get_conductivity_interpolator()
        self.base_phi_interp = self.pore_obj.get_phi_interpolator()
        self.base_dist_interp = self.pore_obj.get_distance_interpolator()
        
        if self.pore_type == "bin_file":
            self.bin_dimensions = self.pore_obj.get_dimensions()
            self.bin_grid_shape = self.pore_obj.get_grid_shape()
        
        if self.rank == 0:
            logger.info("Base conductivity interpolator created")
    
    def _create_cylindrical_conductivity_grid(self, X, Y, Z):
        """Create conductivity grid for cylindrical pore."""
        # Initialize with bulk conductivity
        base_conductivity = np.ones(X.shape) * self.bulk_conductivity
        
        membrane_half_thickness = self.membrane_thickness / 2
        pore_radius = self.pore_radius
        corner_radius = self.corner_radius
        
        membrane_mask = PoreGeometry.create_cylindrical_pore(
            X, Y, Z, pore_radius, membrane_half_thickness, corner_radius
        )
        
        # Set low conductivity in membrane regions
        base_conductivity[membrane_mask] = self.membrane_conductivity
        
        return base_conductivity
    
    def _create_double_cone_conductivity_grid(self, X, Y, Z):
        """Create conductivity grid for double cone pore."""
        # Initialize with bulk conductivity
        base_conductivity = np.ones(X.shape) * self.bulk_conductivity
        
        membrane_half_thickness = self.membrane_thickness / 2
        inner_radius = self.pore_radius
        outer_radius = self.outer_radius
        
        membrane_mask = PoreGeometry.create_double_cone_pore(
            X, Y, Z, inner_radius, outer_radius, membrane_half_thickness
        )
        
        # Set low conductivity in membrane regions
        base_conductivity[membrane_mask] = self.membrane_conductivity
        
        return base_conductivity
    

    def _create_biological_conductivity_grid(self, X, Y, Z):
        """Create conductivity grid for biological pore."""
        membrane_half_thickness = self.membrane_thickness / 2
        
        conductivity_grid, phi_interp = PoreGeometry.create_biological_pore(
                X, Y, Z, 
                pore_pdb=self.biological_pore_pdb,
                membrane_half_thickness=membrane_half_thickness,
                bulk_conductivity=self.bulk_conductivity,
                # cutoff=self.cutoff,
                use_vdw_radii=self.use_vdw_radii,
                default_radius=self.default_radius,
                membrane_conductivity=self.membrane_conductivity,
                membrane_z_offset = self.membrane_z_offset,
                resolution=self.grid_resolution,
                use_charges=self.use_charges,  # NEW
                debye_length=self.debye_length,  # NEW
                bjerrum_length=self.bjerrum_length,  # NEW
                cleanup_temp_files=True,
                box_dimensions=self.box_dimensions,  # Pass the calculated box_dimensions
                use_pdb2pqr = self.use_pdb2pqr
            )
        return conductivity_grid, phi_interp

    def _resolve_gmsh_fine_center(self):
        if self.gmsh_fine_center is not None:
            return
        if self.gmsh_fine_center_mode == "analyte_com":
            if self.prepare_analyte and self.moving_universe is not None:
                self.gmsh_fine_center = [float(v) for v in self.moving_com]
            else:
                if self.rank == 0:
                    logger.warning(
                        "gmsh_fine_center_mode=analyte_com requested, but analyte is unavailable. "
                        "Falling back to origin."
                    )
                self.gmsh_fine_center = [0.0, 0.0, 0.0]
        else:
            self.gmsh_fine_center = [0.0, 0.0, 0.0]

    def _apply_gmsh_options(self, gmsh):
        options = [
            ("General.NumThreads", self.gmsh_num_threads),
            ("Mesh.RandomSeed", self.gmsh_random_seed),
            ("Mesh.RandomFactor", self.gmsh_random_factor),
        ]
        for name, value in options:
            if value is None:
                continue
            try:
                gmsh.option.setNumber(name, float(value))
            except Exception as exc:
                if self.rank == 0 and name not in self._gmsh_option_warned:
                    logger.warning("Gmsh option %s not applied: %s", name, exc)
                self._gmsh_option_warned.add(name)

    def _write_mesh_xdmf(self):
        if not self.save_mesh_xdmf:
            return
        suffix = ""
        if self.update_mesh_each_step:
            suffix = f"_{self._mesh_write_count:04d}"
        path = f"{self.output_prefix}_mesh{suffix}.xdmf"
        try:
            with io.XDMFFile(self.comm, path, "w") as xdmf:
                xdmf.write_mesh(self.mesh)
            if self.rank == 0 and self.verbose_output:
                logger.info("Saved mesh to %s", path)
        except Exception as exc:
            if self.rank == 0:
                logger.warning("Failed to write mesh XDMF (%s): %s", path, exc)
        self._mesh_write_count += 1

    def _maybe_rebuild_mesh_for_position(self, z_position):
        if not self.update_mesh_each_step:
            return 0.0
        if not self.prepare_analyte or self.moving_com is None:
            if self.rank == 0 and not self._gmsh_center_warned:
                logger.warning(
                    "gmsh_fine_center_mode=analyte_com requested, but analyte COM is unavailable. "
                    "Mesh updates will be skipped."
                )
            self._gmsh_center_warned = True
            return 0.0
        new_center = np.array([self.moving_com[0], self.moving_com[1], z_position], dtype=float)
        if self._gmsh_mesh_center is not None and np.allclose(new_center, self._gmsh_mesh_center):
            return 0.0
        self.gmsh_fine_center = [float(v) for v in new_center]
        if self.rank == 0 and self.verbose_output:
            logger.info("Rebuilding Gmsh mesh with fine center at (%.2f, %.2f, %.2f) Å",
                        *self.gmsh_fine_center)
        rebuild_start = time.time()
        self.setup_dolfinx()
        self._gmsh_mesh_center = new_center
        return time.time() - rebuild_start

    def _create_gmsh_mesh(self, extent_x, extent_y, extent_z):
        """Create a Gmsh mesh with optional center refinement box."""
        try:
            import importlib
            import sys
            gmsh = importlib.import_module("gmsh")
            if not hasattr(gmsh, "initialize"):
                sys.modules.pop("gmsh", None)
                gmsh = importlib.import_module("gmsh")
            if not hasattr(gmsh, "initialize"):
                module_path = getattr(gmsh, "__file__", "unknown")
                raise RuntimeError(
                    "Imported module 'gmsh' does not expose initialize "
                    f"(got {gmsh.__name__} from {module_path}). "
                    "Ensure the gmsh package is installed and not shadowed."
                )
        except ImportError as exc:
            raise RuntimeError(
                "mesh_engine='gmsh' requires the gmsh Python package. Install with `pip install gmsh`."
            ) from exc
        try:
            from dolfinx.io import gmsh as dolfinx_gmsh
        except ImportError:
            try:
                import dolfinx.io.gmsh as dolfinx_gmsh
            except ImportError as exc:
                raise RuntimeError(
                    "mesh_engine='gmsh' requires dolfinx gmsh support (gmsh). "
                    "Please update dolfinx or install a build that includes gmsh."
                ) from exc

        self._resolve_gmsh_fine_center()
        fine_size = self.gmsh_fine_size if self.gmsh_fine_size is not None else self.grid_resolution
        coarse_size = self.gmsh_coarse_size if self.gmsh_coarse_size is not None else max(self.grid_resolution * 3.0, self.grid_resolution)

        if self.rank == 0:
            gmsh.initialize()
            self._apply_gmsh_options(gmsh)
            gmsh.model.add("sem_mesh")

            x0 = -extent_x / 2.0
            y0 = -extent_y / 2.0
            z0 = -extent_z / 2.0
            vol = gmsh.model.occ.addBox(x0, y0, z0, extent_x, extent_y, extent_z)
            gmsh.model.occ.synchronize()
            gmsh.model.addPhysicalGroup(3, [vol], 1)
            gmsh.model.setPhysicalName(3, 1, "domain")

            if self.gmsh_fine_box is not None:
                hx, hy, hz = self.gmsh_fine_box
                cx, cy, cz = self.gmsh_fine_center
                gmsh.model.mesh.field.add("Box", 1)
                gmsh.model.mesh.field.setNumber(1, "VIn", fine_size)
                gmsh.model.mesh.field.setNumber(1, "VOut", coarse_size)
                gmsh.model.mesh.field.setNumber(1, "XMin", cx - hx)
                gmsh.model.mesh.field.setNumber(1, "XMax", cx + hx)
                gmsh.model.mesh.field.setNumber(1, "YMin", cy - hy)
                gmsh.model.mesh.field.setNumber(1, "YMax", cy + hy)
                gmsh.model.mesh.field.setNumber(1, "ZMin", cz - hz)
                gmsh.model.mesh.field.setNumber(1, "ZMax", cz + hz)
                gmsh.model.mesh.field.setAsBackgroundMesh(1)
            else:
                gmsh.option.setNumber("Mesh.CharacteristicLengthMin", coarse_size)
                gmsh.option.setNumber("Mesh.CharacteristicLengthMax", coarse_size)

            gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
            gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
            gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
            gmsh.model.mesh.generate(3)

        mesh_data = dolfinx_gmsh.model_to_mesh(gmsh.model, self.comm, 0, gdim=3)
        mesh = mesh_data.mesh

        if self.rank == 0:
            gmsh.finalize()

        return mesh
    
    def setup_dolfinx(self):
        """Initialize DOLFINx solver components."""
        if self.rank == 0:
            logger.info("Creating DOLFINx mesh...")
        
        if self.pore_type == "bin_file":
            # Use original code approach for bin files
            if self.rank == 0:
                logger.info("Using bin file dimensions for mesh sizing")
                if self.xy_margin != 0.0:
                    logger.info("Applying XY margin of %.3f Å to mesh extent", self.xy_margin)

            if self.bin_dimensions is None or self.bin_grid_shape is None:
                raise RuntimeError("Binary pore dimensions are not available; ensure bin file was loaded correctly.")

            base_extent_x, base_extent_y, base_extent_z = self.bin_dimensions
            extent_x = base_extent_x + self.xy_margin
            extent_y = base_extent_y + self.xy_margin
            extent_z = base_extent_z

            if extent_x <= 0 or extent_y <= 0:
                raise ValueError("xy_margin results in non-positive XY extent. Please choose a smaller margin.")

            # Determine cell spacing from the underlying grid
            def _spacing(total_extent: float, cells: int) -> float:
                return total_extent / max(cells - 1, 1)

            spacing_x = _spacing(base_extent_x, self.bin_grid_shape[0])
            spacing_y = _spacing(base_extent_y, self.bin_grid_shape[1])
            spacing_z = _spacing(base_extent_z, self.bin_grid_shape[2])

            nx = max(int(np.ceil(extent_x / spacing_x)), 1)
            ny = max(int(np.ceil(extent_y / spacing_y)), 1)
            nz = max(int(np.ceil(extent_z / spacing_z)), 1)

            if self.mesh_engine == "gmsh":
                self.mesh = self._create_gmsh_mesh(extent_x, extent_y, extent_z)
            else:
                self.mesh = mesh.create_box(
                    self.comm,
                    np.array([[-extent_x / 2., -extent_y / 2., -extent_z / 2.],
                              [extent_x / 2., extent_y / 2., extent_z / 2.]]),
                    [nx, ny, nz]
                )
            
            if self.rank == 0:
                if self.mesh_engine == "gmsh":
                    logger.info("Created Gmsh bin file mesh (reference grid %dx%dx%d)", nx, ny, nz)
                else:
                    logger.info(f"Created bin file mesh with {nx}x{ny}x{nz} elements")
                logger.info(f"Mesh size: {extent_x:.1f} x {extent_y:.1f} x {extent_z:.1f} Å")
            
            sizex, sizey, sizez = extent_x, extent_y, extent_z
        else:
            # Use modern approach for other pore types - MATCH FEniCS EXACTLY
            # Calculate mesh dimensions based on box dimensions and resolution
            base_extent_x = self.box_dimensions['x'][1] - self.box_dimensions['x'][0]
            base_extent_y = self.box_dimensions['y'][1] - self.box_dimensions['y'][0]
            base_extent_z = self.box_dimensions['z'][1] - self.box_dimensions['z'][0]

            extent_x = base_extent_x + self.xy_margin
            extent_y = base_extent_y + self.xy_margin
            extent_z = base_extent_z

            if extent_x <= 0 or extent_y <= 0:
                raise ValueError("xy_margin results in non-positive XY extent. Please choose a smaller margin.")

            nx = max(int(np.ceil(extent_x / self.grid_resolution)), 1)
            ny = max(int(np.ceil(extent_y / self.grid_resolution)), 1)
            nz = max(int(np.ceil(extent_z / self.grid_resolution)), 1)
            
            if self.mesh_engine == "gmsh":
                self.mesh = self._create_gmsh_mesh(extent_x, extent_y, extent_z)
            else:
                # Create mesh exactly like FEniCS version for consistency
                self.mesh = mesh.create_box(
                    self.comm,
                    np.array([[-extent_x / 2., -extent_y / 2., -extent_z / 2.],
                              [extent_x / 2., extent_y / 2., extent_z / 2.]]),
                    [nx, ny, nz]
                )

            if self.rank == 0:
                if self.mesh_engine == "gmsh":
                    logger.info("Created Gmsh mesh (reference grid %dx%dx%d)", nx, ny, nz)
                else:
                    logger.info(f"Created mesh with {nx}x{ny}x{nz} elements")
                logger.info(f"Mesh size: {extent_x:.1f} x {extent_y:.1f} x {extent_z:.1f} Å")

            sizex, sizey, sizez = extent_x, extent_y, extent_z
        
        # Define function spaces
        self.V = fem.functionspace(self.mesh, ("Lagrange", 1))  # For potential
        self.Q = fem.functionspace(self.mesh, ("DG", 0))  # For conductivity
        
        # Define boundary conditions using DOLFINx approach
        fdim = self.mesh.topology.dim - 1
        
        # Calculate z bounds exactly like FEniCS version
        z_min = -sizez / 2.
        z_max = sizez / 2.
        
        # Find boundary facets. Use tolerance tied to grid spacing so coarse meshes
        # still capture the outer planes.
        z_spacing = extent_z / nz if nz > 0 else extent_z
        tol = max(1e-14, 1e-6 * z_spacing)
        top_facets = dolfinx.mesh.locate_entities_boundary(
            self.mesh, fdim, lambda x: np.abs(x[2] - z_max) < tol)
        
        bot_facets = dolfinx.mesh.locate_entities_boundary(
            self.mesh, fdim, lambda x: np.abs(x[2] - z_min) < tol)
        
        # Create boundary conditions
        top_dofs = fem.locate_dofs_topological(self.V, fdim, top_facets)
        bot_dofs = fem.locate_dofs_topological(self.V, fdim, bot_facets)
        
        self.bc_top = fem.dirichletbc(PETSc.ScalarType(self.voltage), top_dofs, self.V)
        self.bc_bot = fem.dirichletbc(PETSc.ScalarType(0.0), bot_dofs, self.V)
        self.bcs = [self.bc_top, self.bc_bot]
        
        # Create facet tags for boundary marking
        facet_indices = np.concatenate([top_facets, bot_facets])
        facet_values = np.concatenate([
            np.full_like(top_facets, 1, dtype=np.int32),  # Top = 1
            np.full_like(bot_facets, 2, dtype=np.int32)   # Bottom = 2
        ])
        
        self.facet_tag = dolfinx.mesh.meshtags(self.mesh, fdim, facet_indices, facet_values)
        
        # Define variational problem
        u = ufl.TrialFunction(self.V)
        v = ufl.TestFunction(self.V)
        self.sig = fem.Function(self.Q)  # Conductivity function
        
        # Define the variational form
        self.a = self.sig * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
        self.L = fem.Constant(self.mesh, PETSc.ScalarType(0.0)) * v * ufl.dx
        
        self.domain_min = [-sizex / 2., -sizey / 2., -sizez / 2.]
        self.domain_max = [sizex / 2., sizey / 2., sizez / 2.]
        self.num_cells = [nx, ny, nz]

        if self.rank == 0:
            logger.info("DOLFINx setup complete")
        self._write_mesh_xdmf()
    
    def get_conductivity_at_position(self, z_position):
        """
        Calculate conductivity field when moving atoms are at given z position.
        Enhanced to use van der Waals radii and support all pore types including binary files.
        Uses the loadFunc approach from original code for better consistency.
        
        Args:
            z_position: Z coordinate to place the center of mass of moving atoms
            
        Returns:
            Conductivity field as numpy array
        """
        # Get base conductivity at mesh points using loadFunc approach
        # moving_atoms = self.moving_universe.atoms
        
        # Calculate current COM of moving atoms
        # moving_com = moving_atoms.center_of_mass()
        displacement = np.array([0, 0, z_position - self.moving_com[2]])
        
        # Get displaced positions in Angstroms after applying rotation
        rotated_positions = self._get_rotated_analyte_positions()
        moving_positions = rotated_positions + displacement
        
        # Get radii in Angstroms (use the pre-assigned VdW radii)
        moving_radii = self.moving_radii

        if self.prevent_analyte_overlap and self.use_radius_overlap_check:
            self._assert_radius_overlap(moving_positions, moving_radii)
        
        # Use loadFunc approach like original code
        # First, load base conductivity
        loadFunc(self.mesh, self.Q, self.sig, self.base_cond_interp, self.bulk_conductivity)
        
        # Get the loaded values
        base_values = self.sig.x.array[:]
        
        # Get mesh coordinates for analyte modification
        mesh_coords = get_dof_coordinates(self.mesh, self.Q)
        
        # Calculate modification due to analyte
        analyte_cond = self.calculate_analyte_conductivity_modification(
            mesh_coords, moving_positions, moving_radii, base_values
        )
        
        # Load modified conductivity back
        self.sig.x.array[:] = analyte_cond
        self.sig.x.scatter_forward()

        if self.prevent_analyte_overlap and self.rank == 0:
            logger.info("Overlap checks passed at Z=%.2f Å.", z_position)

        return analyte_cond

    @staticmethod
    def _distance_to_membrane(R, abs_z, local_radius, membrane_half_thickness):
        radial_term = np.maximum(local_radius - R, 0.0)
        vertical_term = np.maximum(abs_z - membrane_half_thickness, 0.0)
        return np.sqrt(radial_term**2 + vertical_term**2)

    def _assert_radius_overlap(self, atom_positions, atom_radii):
        """Ensure analyte hard-core spheres do not overlap pore solids."""
        if not self.prevent_analyte_overlap or not self.use_radius_overlap_check:
            return
        if atom_positions is None or len(atom_positions) == 0:
            return

        buffer = self.overlap_buffer
        fixed_threshold = self.overlap_distance_threshold
        atom_radii = np.asarray(atom_radii, dtype=float)

        if self.pore_type in ("cylindrical", "double_cone"):
            membrane_half_thickness = self.membrane_thickness / 2.0
            if membrane_half_thickness <= 0:
                return
            R = np.sqrt(atom_positions[:, 0] ** 2 + atom_positions[:, 1] ** 2)
            abs_z = np.abs(atom_positions[:, 2])
            if self.pore_type == "cylindrical":
                local_radius = np.full_like(R, self.pore_radius, dtype=float)
            else:
                z_fraction = np.clip(abs_z / membrane_half_thickness, 0.0, 1.0)
                local_radius = self.pore_radius + (self.outer_radius - self.pore_radius) * z_fraction
            distances = self._distance_to_membrane(R, abs_z, local_radius, membrane_half_thickness)
            if fixed_threshold is not None:
                overlap_mask = distances <= (fixed_threshold + buffer)
            else:
                overlap_mask = distances <= (atom_radii + buffer)
            if np.any(overlap_mask):
                idx = np.flatnonzero(overlap_mask)[0]
                bad_point = atom_positions[idx]
                if fixed_threshold is not None:
                    raise AnalyteOverlapError(
                        "Analyte overlaps membrane wall: "
                        f"distance {distances[idx]:.3f} Å <= threshold "
                        f"{fixed_threshold:.3f} Å + buffer {buffer:.3f} Å at "
                        f"({bad_point[0]:.2f}, {bad_point[1]:.2f}, {bad_point[2]:.2f}) Å."
                    )
                raise AnalyteOverlapError(
                    "Analyte hard core overlaps membrane wall: "
                    f"distance {distances[idx]:.3f} Å <= radius {atom_radii[idx]:.3f} Å "
                    f"+ buffer {buffer:.3f} Å at ({bad_point[0]:.2f}, "
                    f"{bad_point[1]:.2f}, {bad_point[2]:.2f}) Å."
                )
            return

        if self.pore_type == "conical":
            membrane_half_thickness = self.membrane_thickness / 2.0
            if membrane_half_thickness <= 0:
                return
            if self.top_radius is None or self.bottom_radius is None:
                raise AnalyteOverlapError(
                    "Conical overlap check requires both top_radius and bottom_radius "
                    "to be set on the SEM instance."
                )
            R = np.sqrt(atom_positions[:, 0] ** 2 + atom_positions[:, 1] ** 2)
            signed_z = atom_positions[:, 2]
            abs_z = np.abs(signed_z)
            # Asymmetric linear interpolation in *signed* z, matching
            # ConicalPore.get_conductivity_interpolator in pore_geometry.py:
            #   t = 0 at z = -half_thickness (bottom face)
            #   t = 1 at z = +half_thickness (top face)
            thickness = 2.0 * membrane_half_thickness
            t = np.clip((signed_z + membrane_half_thickness) / thickness, 0.0, 1.0)
            local_radius = self.bottom_radius + (self.top_radius - self.bottom_radius) * t
            distances = self._distance_to_membrane(R, abs_z, local_radius, membrane_half_thickness)
            if fixed_threshold is not None:
                overlap_mask = distances <= (fixed_threshold + buffer)
            else:
                overlap_mask = distances <= (atom_radii + buffer)
            if np.any(overlap_mask):
                idx = np.flatnonzero(overlap_mask)[0]
                bad_point = atom_positions[idx]
                if fixed_threshold is not None:
                    raise AnalyteOverlapError(
                        "Analyte overlaps conical membrane wall: "
                        f"distance {distances[idx]:.3f} Å <= threshold "
                        f"{fixed_threshold:.3f} Å + buffer {buffer:.3f} Å at "
                        f"({bad_point[0]:.2f}, {bad_point[1]:.2f}, {bad_point[2]:.2f}) Å "
                        f"(local pore radius {local_radius[idx]:.3f} Å)."
                    )
                raise AnalyteOverlapError(
                    "Analyte hard core overlaps conical membrane wall: "
                    f"distance {distances[idx]:.3f} Å <= radius {atom_radii[idx]:.3f} Å "
                    f"+ buffer {buffer:.3f} Å at ({bad_point[0]:.2f}, "
                    f"{bad_point[1]:.2f}, {bad_point[2]:.2f}) Å "
                    f"(local pore radius {local_radius[idx]:.3f} Å)."
                )
            return

        if self.pore_type == "biological":
            pore_positions = getattr(self.pore_obj, "pore_positions", None)
            pore_radii = getattr(self.pore_obj, "pore_radii", None)
            pore_tree = getattr(self.pore_obj, "pore_tree", None)
            if pore_positions is None or pore_radii is None or len(pore_positions) == 0:
                return
            if pore_tree is None:
                pore_tree = KDTree(pore_positions)
                self.pore_obj.pore_tree = pore_tree
            distances, indices = pore_tree.query(atom_positions)
            if fixed_threshold is not None:
                overlap_mask = distances <= (fixed_threshold + buffer)
            else:
                overlap_mask = distances <= (atom_radii + pore_radii[indices] + buffer)
            if np.any(overlap_mask):
                idx = np.flatnonzero(overlap_mask)[0]
                bad_point = atom_positions[idx]
                if fixed_threshold is not None:
                    raise AnalyteOverlapError(
                        "Analyte overlaps biological pore atoms: "
                        f"distance {distances[idx]:.3f} Å <= threshold "
                        f"{fixed_threshold:.3f} Å + buffer {buffer:.3f} Å at "
                        f"({bad_point[0]:.2f}, {bad_point[1]:.2f}, {bad_point[2]:.2f}) Å."
                    )
                raise AnalyteOverlapError(
                    "Analyte hard core overlaps biological pore atoms: "
                    f"distance {distances[idx]:.3f} Å <= radius sum "
                    f"{atom_radii[idx] + pore_radii[indices[idx]]:.3f} Å + buffer {buffer:.3f} Å "
                    f"at ({bad_point[0]:.2f}, {bad_point[1]:.2f}, {bad_point[2]:.2f}) Å."
                )
            return

        if self.pore_type == "bin_file":
            if self.base_dist_interp is None:
                if self.rank == 0:
                    logger.warning(
                        "Radius-sum overlap check skipped: bin_file has no distance map."
                    )
                return
            try:
                dist_vals = np.asarray(self.base_dist_interp(atom_positions))
            except ValueError as exc:
                raise AnalyteOverlapError(
                    f"Analyte positions outside distance interpolation domain: {exc}"
                ) from exc
            if dist_vals.ndim == 0:
                dist_vals = dist_vals.reshape(1)
            if fixed_threshold is not None:
                overlap_mask = ~np.isfinite(dist_vals) | (dist_vals <= (fixed_threshold + buffer))
            else:
                overlap_mask = ~np.isfinite(dist_vals) | (dist_vals <= (atom_radii + buffer))
            if np.any(overlap_mask):
                idx = np.flatnonzero(overlap_mask)[0]
                bad_point = atom_positions[idx]
                if fixed_threshold is not None:
                    raise AnalyteOverlapError(
                        "Analyte overlaps bin_file pore region: "
                        f"distance {dist_vals[idx]:.3f} Å <= threshold "
                        f"{fixed_threshold:.3f} Å + buffer {buffer:.3f} Å at "
                        f"({bad_point[0]:.2f}, {bad_point[1]:.2f}, {bad_point[2]:.2f}) Å."
                    )
                raise AnalyteOverlapError(
                    "Analyte hard core overlaps bin_file pore region: "
                    f"distance {dist_vals[idx]:.3f} Å <= radius {atom_radii[idx]:.3f} Å "
                    f"+ buffer {buffer:.3f} Å at ({bad_point[0]:.2f}, "
                    f"{bad_point[1]:.2f}, {bad_point[2]:.2f}) Å."
                )

    def _get_rotated_analyte_positions(self):
        """Return analyte coordinates with the current rotation applied."""
        if self._base_moving_positions is None or len(self._base_moving_positions) == 0:
            return self.moving_positions
        centered = self._base_moving_positions - self.moving_com
        rotated = centered @ self._current_rotation_matrix.T + self.moving_com
        return rotated

    def set_analyte_rotation_matrix(self, rotation_matrix: np.ndarray | None):
        """Set the analyte rotation matrix (3x3)."""
        if rotation_matrix is None:
            self._current_rotation_matrix = np.eye(3)
            return
        matrix = np.asarray(rotation_matrix, dtype=float)
        if matrix.shape != (3, 3):
            raise ValueError("rotation_matrix must be 3x3")
        self._current_rotation_matrix = matrix

    def reset_analyte_rotation(self):
        """Reset analyte to its original orientation."""
        self._current_rotation_matrix = np.eye(3)
    
    def calculate_analyte_conductivity_modification(self, mesh_coords, atom_positions, 
                                                   atom_radii, base_conductivity):
        """
        Modify conductivity based on presence of analyte atoms.
        Now uses accurate van der Waals radii for each atom.
        
        Args:
            mesh_coords: Nx3 array of mesh coordinates (Å)
            atom_positions: Mx3 array of atom positions (Å)
            atom_radii: M array of atom radii (Å) - now element-specific VdW radii
            base_conductivity: N array of base conductivity values
            
        Returns:
            N array of modified conductivity values (S/m)
        """
        # Start with base conductivity
        conductivity = base_conductivity.copy()
        
        # Build KDTree for efficient distance queries
        atom_tree = KDTree(atom_positions)
        
        # Query distances
        cutoff = self.cutoff  # Keep in Angstroms
        distances, indices = atom_tree.query(
            mesh_coords, 
            distance_upper_bound=cutoff
        )
        
        # Adjust distances for atom radii
        valid_mask = ~np.isinf(distances)
        if np.any(valid_mask):
            distances[valid_mask] -= atom_radii[indices[valid_mask]]
            distances[valid_mask] = np.maximum(distances[valid_mask], 0)

            # Apply conductivity model to get modulation factor
            analyte_cond = self.conductivity_model(distances[valid_mask])

            if self.prevent_analyte_overlap:
                min_distance = getattr(self.conductivity_model, "min_distance", 0.0)
                hard_core_mask = distances[valid_mask] <= min_distance
                if np.any(hard_core_mask):
                    base_vals = base_conductivity[valid_mask]
                    membrane_mask = base_vals <= self.membrane_conductivity
                    overlap_mask = hard_core_mask & membrane_mask
                    if np.any(overlap_mask):
                        valid_idx = np.flatnonzero(valid_mask)
                        first_local = np.flatnonzero(overlap_mask)[0]
                        global_idx = valid_idx[first_local]
                        bad_point = mesh_coords[global_idx]
                        bad_base = float(base_conductivity[global_idx])
                        bad_analyte = float(analyte_cond[first_local])
                        raise AnalyteOverlapError(
                            "Analyte hard core overlaps membrane/pore region: "
                            f"{overlap_mask.sum()} mesh points (base {bad_base:.3e} S/m, "
                            f"analyte {bad_analyte:.3e} S/m). First at "
                            f"({bad_point[0]:.2f}, {bad_point[1]:.2f}, {bad_point[2]:.2f}) Å."
                        )

            # Modify conductivity: use minimum of base and analyte-modified conductivity
            conductivity[valid_mask] = np.minimum(
                conductivity[valid_mask],
                analyte_cond
            )
        
        return conductivity
    
    def solve_for_current(self, conductivity=None):
        """
        Solve FEM problem for given conductivity field and return current.
        Modified to work with DOLFINx.
        
        Args:
            conductivity: Conductivity field values at conductivity DOFs (optional, uses current sig if None)
            
        Returns:
            current: Calculated current (A)
        """
        try:
            if self.rank == 0:
                logger.info("Starting solve_for_current (DOLFINx)...")
            
            # If conductivity provided, load it
            if conductivity is not None:
                if self.rank == 0:
                    logger.info("Loading provided conductivity...")
                self.sig.x.array[:] = conductivity
                self.sig.x.scatter_forward()
                if self.rank == 0:
                    logger.info("Conductivity loaded")
            else:
                if self.rank == 0:
                    logger.info("Using existing conductivity in self.sig")
            
            # Solve using LinearProblem
            if self.rank == 0:
                logger.info("Setting up and solving linear problem...")
            
            problem = LinearProblem(
                self.a, self.L, bcs=self.bcs,
                petsc_options={
                    "ksp_type": "gmres",
                    "pc_type": "hypre",
                    "ksp_rtol": 1e-10,
                    "ksp_max_it": 40000
                },
                petsc_options_prefix="sem_"
            )
            
            uh = problem.solve()
            
            if self.rank == 0:
                logger.info("System solved")
            
            # Calculate flux at boundaries - MATCH FEniCS EXACTLY
            if self.rank == 0:
                logger.info("Calculating flux at boundaries...")
            
            # Use fixed normal vector like FEniCS version for consistency
            fixed_normal = ufl.as_vector([0.0, 0.0, 1.0])  # Match FEniCS Constant((0,0,1))
            ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.facet_tag)
            
            flux_form_top = ufl.dot(fixed_normal, self.sig * ufl.grad(uh)) * ds(1)
            flux_form_bot = ufl.dot(fixed_normal, self.sig * ufl.grad(uh)) * ds(2)
            
            flux_top = fem.assemble_scalar(fem.form(flux_form_top))
            flux_bot = fem.assemble_scalar(fem.form(flux_form_bot))
            
            # Gather results across MPI processes
            if self.comm.size > 1:
                flux_top = self.comm.allreduce(flux_top, op=MPI.SUM)
                flux_bot = self.comm.allreduce(flux_bot, op=MPI.SUM)
            
            if self.rank == 0:
                logger.info(f"Flux top: {flux_top}, Flux bottom: {flux_bot}")
            
            # Return top flux like original
            current = abs(flux_top)
            if self.rank == 0:
                logger.info(f"Returning current: {current}")
            return current
            
        except Exception as e:
            if self.rank == 0:
                logger.error(f"Error in solve_for_current: {e}")
                logger.error(f"Error type: {type(e)}")
                import traceback
                logger.error(f"Traceback: {traceback.format_exc()}")
            raise
    
    def calculate_open_pore_current(self):
        """Calculate baseline current with no analyte in the system."""
        if self.rank == 0:
            logger.info("Calculating open pore current...")
        
        try:
            # Load base conductivity using loadFunc
            if self.rank == 0:
                logger.info("Loading base conductivity into DOLFINx function...")
            loadFunc(self.mesh, self.Q, self.sig, self.base_cond_interp, self.bulk_conductivity)
            if self.rank == 0:
                logger.info("Base conductivity loaded successfully")
            
            # Solve for current
            if self.rank == 0:
                logger.info("Solving for current...")
            open_current = self.solve_for_current()
            if self.rank == 0:
                logger.info(f"Open pore current: {open_current:.6e} nA")
            self._open_pore_current = open_current
            return open_current
            
        except Exception as e:
            if self.rank == 0:
                logger.error(f"Error in calculate_open_pore_current: {e}")
                logger.error(f"Error type: {type(e)}")
                import traceback
                logger.error(f"Traceback: {traceback.format_exc()}")
            raise
    
    def run(self, open_current=None):
        """
        Run the complete vertical movement simulation with robust timing.
        
        Args:
            open_current: Optional precomputed open-pore current to reuse across runs.

        Returns:
            results: Dictionary with z_positions, currents, normalized currents, and timing data
        """
        # Start total simulation timer
        total_start_time = time.time()

        # If gmsh fine center follows analyte, align mesh to first position before open-pore solve
        initial_mesh_time = self._maybe_rebuild_mesh_for_position(self.z_start)
        
        # First calculate open pore current
        open_pore_time = 0.0
        reused_open_pore = False
        if open_current is None:
            if self._open_pore_current is not None:
                open_current = self._open_pore_current
                reused_open_pore = True
            else:
                if self.rank == 0:
                    logger.info("Calculating open pore current...")
                open_pore_start = time.time()
                open_current = self.calculate_open_pore_current()
                open_pore_time = time.time() - open_pore_start
        else:
            open_current = float(open_current)
            self._open_pore_current = open_current
            reused_open_pore = True
        
        if self.rank == 0:
            if reused_open_pore:
                logger.info(f"Using cached open pore current: {open_current:.6e} nA")
            else:
                logger.info(f"Open pore current calculated in {open_pore_time:.2f} seconds: {open_current:.6e} nA")
        
        # Calculate number of steps
        num_steps = int(abs(self.z_end - self.z_start) / self.z_step) + 1
        direction = np.sign(self.z_end - self.z_start)
        z_positions = self.z_start + np.arange(num_steps) * (direction * self.z_step)
        
        if self.rank == 0:
            logger.info(f"Running simulation with {num_steps} positions")
            logger.info(f"Z range: {self.z_start} to {self.z_end} Å, step: {self.z_step} Å")
            logger.info(f"Pore type: {self.pore_type}")
            if self.pore_type == "cylindrical":
                logger.info(f"Pore radius: {self.pore_radius} Å")
                if self.corner_radius > 0:
                    logger.info(f"Corner radius: {self.corner_radius} Å")
            elif self.pore_type == "double_cone":
                logger.info(f"Inner radius: {self.pore_radius} Å")
                logger.info(f"Outer radius: {self.outer_radius} Å")
            elif self.pore_type == "biological":
                logger.info(f"Biological pore: {self.biological_pore_pdb}")
            elif self.pore_type == "bin_file":
                logger.info(f"Binary file: {self.bin_file_path}")
            
            if self.use_vdw_radii:
                logger.info("Using van der Waals radii for accurate atomic volumes")
            else:
                logger.info(f"Using uniform radius of {self.default_radius} Å for all atoms")
        
        # Initialize arrays for results and timing
        currents = []
        position_times = []
        conductivity_times = []
        solver_times = []
        mesh_times = []
        
        # Start main simulation loop
        simulation_start_time = time.time()
        
        for i, z_pos in enumerate(z_positions):
            position_start_time = time.time()
            
            if self.rank == 0 and self.verbose_output:
                logger.info(f"Processing position {i+1}/{num_steps}: Z = {z_pos:.1f} Å")
            
            # Update mesh if fine center follows analyte COM
            mesh_time = self._maybe_rebuild_mesh_for_position(z_pos)
            mesh_times.append(mesh_time)

            try:

                # Time conductivity calculation
                conductivity_start = time.time()
                self.get_conductivity_at_position(z_pos)
                conductivity_time = time.time() - conductivity_start
                conductivity_times.append(conductivity_time)
                
                # Time solver
                solver_start = time.time()
                current = self.solve_for_current()
                solver_time = time.time() - solver_start
                solver_times.append(solver_time)
                
                currents.append(current)
            except AnalyteOverlapError as overlap_exc:
                if self.rank == 0:
                    logger.warning(
                        "Skipping position %d/%d (Z=%.2f Å) due to analyte overlap: %s",
                        i + 1,
                        num_steps,
                        z_pos,
                        overlap_exc,
                    )
                conductivity_times.append(np.nan)
                solver_times.append(np.nan)
                currents.append(np.nan)
                position_times.append(np.nan)
                continue
            
            # Calculate timing for this position
            position_time = time.time() - position_start_time
            position_times.append(position_time)
            
            # Calculate blockage
            blockage = (1 - current/open_current) * 100 if np.isfinite(current) else np.nan
            
            if self.rank == 0 and self.verbose_output:
                logger.info(f"  Current: {current:.6e} nA (blockage: {blockage:.1f}%)")
                logger.info(f"  Timing - Mesh: {mesh_time:.3f}s, "
                        f"Conductivity: {conductivity_time:.3f}s, "
                        f"Solver: {solver_time:.3f}s, Total: {position_time:.3f}s")
                
                # Estimate remaining time (after first few positions for better accuracy)
                if i >= 2:  # Need at least 3 positions for good estimate
                    valid_times = np.array(position_times)[np.isfinite(position_times)]
                    if valid_times.size >= 1:
                        avg_time_per_position = np.mean(valid_times)
                        remaining_positions = num_steps - (i + 1)
                        estimated_remaining_time = avg_time_per_position * remaining_positions
                        
                        hours = int(estimated_remaining_time // 3600)
                        minutes = int((estimated_remaining_time % 3600) // 60)
                        seconds = int(estimated_remaining_time % 60)
                        
                        if hours > 0:
                            time_str = f"{hours}h {minutes}m {seconds}s"
                        elif minutes > 0:
                            time_str = f"{minutes}m {seconds}s"
                        else:
                            time_str = f"{seconds}s"
                        
                        logger.info(f"  Progress: {(i+1)/num_steps*100:.1f}%, "
                                f"Est. remaining: {time_str}")
                    else:
                        logger.info(
                            f"  Progress: {(i+1)/num_steps*100:.1f}%, Est. remaining: n/a"
                        )
            
            # Save intermediate result with timing (only rank 0)
            if self.rank == 0 and self.verbose_output:
                with open(f"{self.output_prefix}_position_{i:04d}.dat", 'w') as f:
                    f.write(f"# Position {i+1}/{num_steps}\n")
                    f.write(f"# Z_position: {z_pos} Å\n")
                    f.write(f"# Current: {current if np.isfinite(current) else float('nan'):.6e} nA\n")
                    f.write(f"# Blockage: {blockage if np.isfinite(blockage) else float('nan'):.2f}%\n")
                    f.write(f"# Mesh_time: {mesh_time:.3f} s\n")
                    f.write(f"# Conductivity_time: {conductivity_time:.3f} s\n")
                    f.write(f"# Solver_time: {solver_time:.3f} s\n")
                    f.write(f"# Total_time: {position_time:.3f} s\n")
                    f.write(f"# Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"{z_pos:.1f} {current:.6e} {blockage:.2f} {mesh_time:.3f} {conductivity_time:.3f} {solver_time:.3f} {position_time:.3f}\n")
        
        # Calculate final timing statistics
        simulation_time = time.time() - simulation_start_time
        total_time = time.time() - total_start_time
        
        # Convert to numpy arrays for statistics
        position_times = np.array(position_times)
        mesh_times = np.array(mesh_times)
        conductivity_times = np.array(conductivity_times)
        solver_times = np.array(solver_times)
        currents = np.array(currents)
        
        normalized_currents = currents / open_current
        blockages = (1 - normalized_currents) * 100
        
        # Calculate comprehensive timing statistics
        if self.rank == 0:
            avg_position_time = np.mean(position_times)
            std_position_time = np.std(position_times)
            avg_mesh_time = np.mean(mesh_times)
            std_mesh_time = np.std(mesh_times)
            avg_conductivity_time = np.mean(conductivity_times)
            std_conductivity_time = np.std(conductivity_times)
            avg_solver_time = np.mean(solver_times)
            std_solver_time = np.std(solver_times)
            
            # Performance metrics
            positions_per_hour = 3600 / avg_position_time if avg_position_time > 0 else 0
            total_calc_time = np.sum(mesh_times) + np.sum(conductivity_times) + np.sum(solver_times)
            overhead_time = simulation_time - total_calc_time
            
            logger.info("=" * 80)
            logger.info("COMPREHENSIVE TIMING ANALYSIS")
            logger.info("=" * 80)
            logger.info(f"Total simulation time: {total_time:.2f} seconds ({total_time/60:.1f} minutes)")
            logger.info(f"Open pore calculation: {open_pore_time:.2f} seconds ({open_pore_time/total_time*100:.1f}%)")
            logger.info(f"Initial mesh rebuild:  {initial_mesh_time:.2f} seconds ({initial_mesh_time/total_time*100:.1f}%)")
            logger.info(f"Main simulation loop: {simulation_time:.2f} seconds ({simulation_time/total_time*100:.1f}%)")
            logger.info(f"Overhead (I/O, logging): {overhead_time:.2f} seconds ({overhead_time/total_time*100:.1f}%)")
            logger.info("")
            logger.info("Per-position timing statistics:")
            logger.info(f"  Mesh rebuild:          {avg_mesh_time:.3f} ± {std_mesh_time:.3f} s")
            logger.info(f"  Conductivity calculation: {avg_conductivity_time:.3f} ± {std_conductivity_time:.3f} s")
            logger.info(f"  FEM solver:              {avg_solver_time:.3f} ± {std_solver_time:.3f} s")
            logger.info(f"  Total per position:      {avg_position_time:.3f} ± {std_position_time:.3f} s")
            logger.info("")
            logger.info("Performance breakdown:")
            logger.info(f"  Conductivity vs Solver:  {avg_conductivity_time/avg_solver_time:.2f}:1 ratio")
            logger.info(f"  Fastest position:        {np.min(position_times):.3f} s")
            logger.info(f"  Slowest position:        {np.max(position_times):.3f} s")
            logger.info(f"  Throughput:              {positions_per_hour:.1f} positions/hour")
            logger.info("")
            logger.info("Time distribution:")
            mesh_pct = np.sum(mesh_times) / total_calc_time * 100
            cond_pct = np.sum(conductivity_times) / total_calc_time * 100
            solver_pct = np.sum(solver_times) / total_calc_time * 100
            logger.info(f"  Mesh rebuilds:           {mesh_pct:.1f}% of compute time")
            logger.info(f"  Conductivity calculations: {cond_pct:.1f}% of compute time")
            logger.info(f"  FEM solver:               {solver_pct:.1f}% of compute time")
            logger.info("=" * 80)
        
        # Save final results with comprehensive timing (only rank 0)
        if self.rank == 0:
            results = {
                'z_positions': z_positions,
                'currents': currents,
                'normalized_currents': normalized_currents,
                'blockages': blockages,
                'open_current': open_current,
                'pore_type': self.pore_type,
                'pore_radius': self.pore_radius,
                'outer_radius': self.outer_radius if self.pore_type == "double_cone" else None,
                'corner_radius': self.corner_radius if self.pore_type == "cylindrical" else None,
                'biological_pore_pdb': self.biological_pore_pdb if self.pore_type == "biological" else None,
                'bin_file_path': self.bin_file_path if self.pore_type == "bin_file" else None,
                'use_vdw_radii': self.use_vdw_radii,
                # Comprehensive timing data
                'timing': {
                    'total_time': total_time,
                    'open_pore_time': open_pore_time,
                    'initial_mesh_time': initial_mesh_time,
                    'simulation_time': simulation_time,
                    'overhead_time': overhead_time,
                    'position_times': position_times,
                    'mesh_times': mesh_times,
                    'conductivity_times': conductivity_times,
                    'solver_times': solver_times,
                    'avg_position_time': avg_position_time,
                    'avg_mesh_time': avg_mesh_time,
                    'avg_conductivity_time': avg_conductivity_time,
                    'avg_solver_time': avg_solver_time,
                    'std_position_time': std_position_time,
                    'std_mesh_time': std_mesh_time,
                    'std_conductivity_time': std_conductivity_time,
                    'std_solver_time': std_solver_time,
                    'positions_per_hour': positions_per_hour,
                    'mesh_percentage': mesh_pct,
                    'conductivity_percentage': cond_pct,
                    'solver_percentage': solver_pct
                }
            }
            
            # Format time for file headers
            def format_time_str(seconds):
                if seconds < 60:
                    return f"{seconds:.1f}s"
                elif seconds < 3600:
                    return f"{seconds/60:.1f}m"
                else:
                    return f"{seconds/3600:.1f}h"
            
            # Save main results file with comprehensive timing columns
            header_lines = [
                f"Z_position(Å) Current(nA) Normalized_Current Blockage(%) Mesh_time(s) Conductivity_time(s) Solver_time(s) Total_time(s)",
                f"Pore_type: {self.pore_type}",
                f"Open_pore_current: {open_current:.6e} A",
                f"Total_simulation_time: {format_time_str(total_time)}",
                f"Average_time_per_position: {avg_position_time:.3f}s",
                f"Throughput: {positions_per_hour:.1f} positions/hour",
                f"Use_VdW_radii: {self.use_vdw_radii}"
            ]
            header = "\n".join([f"# {line}" for line in header_lines])
            
            np.savetxt(f"{self.output_prefix}_results.txt", 
                    np.column_stack([z_positions, currents, normalized_currents, blockages, 
                                    mesh_times, conductivity_times, solver_times, position_times]),
                    header=header, fmt=['%.1f', '%.6e', '%.6f', '%.2f', '%.3f', '%.3f', '%.3f', '%.3f'])
            
            # Save detailed timing analysis
            timing_header_lines = [
                "COMPREHENSIVE TIMING ANALYSIS",
                f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
                f"Total positions: {num_steps}",
                f"Total simulation time: {format_time_str(total_time)}",
                f"Open pore calculation: {format_time_str(open_pore_time)}",
                f"Initial mesh rebuild: {format_time_str(initial_mesh_time)}",
                f"Main simulation: {format_time_str(simulation_time)}",
                f"Average per position: {avg_position_time:.3f}s ± {std_position_time:.3f}s",
                f"Mesh rebuild avg: {avg_mesh_time:.3f}s ± {std_mesh_time:.3f}s",
                f"Conductivity avg: {avg_conductivity_time:.3f}s ± {std_conductivity_time:.3f}s",
                f"Solver avg: {avg_solver_time:.3f}s ± {std_solver_time:.3f}s",
                f"Throughput: {positions_per_hour:.1f} positions/hour",
                "",
                "Position Z_position(Å) Mesh_time(s) Conductivity_time(s) Solver_time(s) Total_time(s)"
            ]
            if self.verbose_output:
                timing_header = "\n".join([f"# {line}" for line in timing_header_lines])
                
                timing_data = np.column_stack([
                    np.arange(1, len(z_positions) + 1),  # Position number
                    z_positions,
                    mesh_times,
                    conductivity_times,
                    solver_times,
                    position_times
                ])
                
                np.savetxt(f"{self.output_prefix}_timing_analysis.txt", timing_data, 
                        header=timing_header, fmt=['%d', '%.1f', '%.6f', '%.6f', '%.6f', '%.6f'])
                
                # Save timing summary for quick reference
                with open(f"{self.output_prefix}_timing_summary.txt", 'w') as f:
                    f.write(f"SEM Simulation Timing Summary\n")
                    f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"{'='*50}\n\n")
                    f.write(f"Total Time: {format_time_str(total_time)}\n")
                    f.write(f"Positions: {num_steps}\n")
                    f.write(f"Throughput: {positions_per_hour:.1f} positions/hour\n")
                    f.write(f"Open Pore: {format_time_str(open_pore_time)} ({open_pore_time/total_time*100:.1f}%)\n")
                    f.write(f"Initial Mesh: {format_time_str(initial_mesh_time)} ({initial_mesh_time/total_time*100:.1f}%)\n")
                    f.write(f"Simulation: {format_time_str(simulation_time)} ({simulation_time/total_time*100:.1f}%)\n\n")
                    f.write(f"Per Position (average):\n")
                    f.write(f"  Mesh:         {avg_mesh_time:.3f}s ({mesh_pct:.1f}%)\n")
                    f.write(f"  Conductivity: {avg_conductivity_time:.3f}s ({cond_pct:.1f}%)\n")
                    f.write(f"  FEM Solver:   {avg_solver_time:.3f}s ({solver_pct:.1f}%)\n")
                    f.write(f"  Total:        {avg_position_time:.3f}s\n\n")
                    f.write(f"Performance:\n")
                    f.write(f"  Fastest: {np.min(position_times):.3f}s\n")
                    f.write(f"  Slowest: {np.max(position_times):.3f}s\n")
                    f.write(f"  Std Dev: {std_position_time:.3f}s\n")
            
            logger.info("Simulation complete!")
            logger.info(f"Maximum blockage: {np.max(blockages):.1f}% at Z = {z_positions[np.argmax(blockages)]:.1f} Å")
            logger.info(f"Results saved to:")
            logger.info(f"  Main results: {self.output_prefix}_results.txt")
            logger.info(f"  Timing analysis: {self.output_prefix}_timing_analysis.txt")
            logger.info(f"  Timing summary: {self.output_prefix}_timing_summary.txt")
            
            return results
        else:
            return None  # Non-root processes don't return results
        
    def get_conductivity_grid_for_preview(self, z_position):
        """
        Get 2D conductivity slice for visualization at given z position.
        Now unified for all pore geometries since they all use the same RegularGridInterpolator.
        
        Args:
            z_position: Z coordinate to place the center of mass of moving atoms
            
        Returns:
            x_coords, z_coords, conductivity_2d: Arrays for plotting
        """
        # Use the exact domain from the DOLFINx mesh
        x_coords = np.linspace(self.domain_min[0], self.domain_max[0], self.num_cells[0] + 1)
        z_coords = np.linspace(self.domain_min[2], self.domain_max[2], self.num_cells[2] + 1)
        y_middle = (self.domain_min[1] + self.domain_max[1]) / 2
        
        # Create meshgrid for the slice
        X_slice, Z_slice = np.meshgrid(x_coords, z_coords, indexing='ij')
        Y_slice = np.full_like(X_slice, y_middle)
        
        # Create coordinate array for interpolation (in Angstroms)
        slice_coords = np.column_stack([
            X_slice.ravel(),  # Keep in Angstroms
            Y_slice.ravel(),  # Keep in Angstroms
            Z_slice.ravel()   # Keep in Angstroms
        ])
        
        # Get base conductivity at slice points using unified interpolator
        base_cond = self.base_cond_interp(slice_coords)
        
        # Handle NaN values
        nan_mask = np.isnan(base_cond)
        base_cond[nan_mask] = self.bulk_conductivity
        
        # Now modify conductivity based on analyte position
        displacement = np.array([0, 0, z_position - self.moving_com[2]])
        
        rotated_positions = self._get_rotated_analyte_positions()
        moving_positions = rotated_positions + displacement
        
        # Get radii in Angstroms (use the pre-assigned VdW radii)
        moving_radii = self.moving_radii
        
        # Calculate modification due to analyte
        analyte_cond = self.calculate_analyte_conductivity_modification(
            slice_coords, moving_positions, moving_radii, base_cond
        )
        
        # Reshape back to 2D
        conductivity_2d = analyte_cond.reshape(X_slice.shape)
        
        return x_coords, z_coords, conductivity_2d
