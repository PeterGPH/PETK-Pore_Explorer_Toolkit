"""
Pore geometry classes and functions.
Handles different pore types: cylindrical, double_cone, biological, and bin_file.
"""

import numpy as np
import logging
import sys
import os
import subprocess
import time
from pathlib import Path
from typing import Optional
from scipy.spatial import KDTree
from scipy.interpolate import RegularGridInterpolator
from abc import ABC, abstractmethod

from .utils import readbinGrid, condfrac
from .van_der_waals import VanDerWaalsRadii
from .conductivity_models import SimpleConductivityModel
from .structure_preparation import prepare_structure, PreparedStructure

logger = logging.getLogger(__name__)

def _conductivity_from_distance(distance_map, bulk_conductivity, membrane_conductivity):
    """
    Convert a gen_dist-style distance field into conductivity.
    Mimics the bin-file workflow: map distance→fraction via condfrac and
    then interpolate between membrane and bulk conductivity.
    """
    fraction = condfrac(distance_map)
    return membrane_conductivity + fraction * (bulk_conductivity - membrane_conductivity)

def _distance_to_membrane(R, abs_Z, local_radius, membrane_half_thickness):
    """
    Analytic distance to the solid membrane defined by:
      R >= local_radius  and  |Z| <= membrane_half_thickness.
    Outside this region we compute Euclidean distance to the allowed set,
    matching the notion of the distance transform produced by gen_dist.
    """
    radial_term = np.maximum(local_radius - R, 0.0)
    vertical_term = np.maximum(abs_Z - membrane_half_thickness, 0.0)
    return np.sqrt(radial_term**2 + vertical_term**2)

class BasePore(ABC):
    """
    Abstract base class for all pore types. Subclasses must implement the interface.
    """
    
    @abstractmethod
    def get_conductivity_interpolator(self):
        """Return a RegularGridInterpolator for conductivity."""
        pass
    
    def get_phi_interpolator(self):
        """Return potential interpolator if available, else a zero-function."""
        return lambda pts: np.zeros(len(pts))  # Default for non-charge-aware pores

    def get_distance_interpolator(self):
        """Return a distance interpolator if available, else None."""
        return None
    
    def get_dimensions(self):
        """Return [Lm, Wm, Hm] or None if not applicable."""
        return None
    
    def get_grid_shape(self):
        """Return [nx, ny, nz] or None if not applicable."""
        return None
    
    def get_grid_spacing(self):
        """Return representative grid spacing in Å if defined."""
        return None

class CylindricalPore(BasePore):
    def __init__(self, X, Y, Z, pore_radius, membrane_half_thickness, 
                 corner_radius=None, chamfer_depth=None, 
                 bulk_conductivity=10.5, membrane_conductivity=0.0001):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.pore_radius = pore_radius
        self.membrane_half_thickness = membrane_half_thickness
        self.corner_radius = corner_radius
        self.chamfer_depth = chamfer_depth
        self.bulk_conductivity = bulk_conductivity
        self.membrane_conductivity = membrane_conductivity
        self._interpolator = None  # Lazy
        self._local_pore_radius = None

    def _compute_local_pore_radius(self):
        base_radius = np.full_like(self.Z, self.pore_radius, dtype=float)
        if self.corner_radius is None or self.corner_radius <= 0:
            return base_radius

        chamfer_depth = self.chamfer_depth if self.chamfer_depth is not None else self.corner_radius
        if chamfer_depth is None or chamfer_depth <= 0:
            return base_radius

        edge_radius = self.pore_radius + self.corner_radius
        z_edge_dist = np.maximum(self.membrane_half_thickness - np.abs(self.Z), 0.0)
        in_chamfer_zone = z_edge_dist < chamfer_depth
        chamfer_progress = np.zeros_like(self.Z, dtype=float)
        chamfer_progress[in_chamfer_zone] = np.clip(
            z_edge_dist[in_chamfer_zone] / chamfer_depth, 0.0, 1.0
        )

        local_radius = np.where(
            in_chamfer_zone,
            edge_radius + (self.pore_radius - edge_radius) * chamfer_progress,
            base_radius,
        )
        return local_radius
    
    def get_conductivity_interpolator(self):
        if self._interpolator is None:
            R = np.sqrt(self.X**2 + self.Y**2)
            local_pore_radius = self._compute_local_pore_radius()
            self._local_pore_radius = local_pore_radius

            distance_map = _distance_to_membrane(
                R,
                np.abs(self.Z),
                local_pore_radius,
                self.membrane_half_thickness,
            )
            conductivity_grid = _conductivity_from_distance(
                distance_map,
                self.bulk_conductivity,
                self.membrane_conductivity,
            )
            
            # Extract edges
            x_range = np.unique(self.X[:, 0, 0])
            y_range = np.unique(self.Y[0, :, 0])
            z_range = np.unique(self.Z[0, 0, :])
            
            self._interpolator = RegularGridInterpolator(
                (x_range, y_range, z_range), conductivity_grid,
                bounds_error=False, fill_value=self.bulk_conductivity
            )
        return self._interpolator

class DoubleConePore(BasePore):
    def __init__(self, X, Y, Z, inner_radius, outer_radius, membrane_half_thickness, 
                 bulk_conductivity=10.5, membrane_conductivity=0.0001):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.membrane_half_thickness = membrane_half_thickness
        self.bulk_conductivity = bulk_conductivity
        self.membrane_conductivity = membrane_conductivity
        self._interpolator = None  # Lazy
    
    def get_conductivity_interpolator(self):
        if self._interpolator is None:
            R = np.sqrt(self.X**2 + self.Y**2)
            abs_z = np.abs(self.Z)
            if self.membrane_half_thickness <= 0:
                raise ValueError("Double-cone pores require a non-zero membrane thickness.")

            z_fraction = np.clip(
                abs_z / self.membrane_half_thickness,
                0.0,
                1.0,
            )
            local_pore_radius = self.inner_radius + (self.outer_radius - self.inner_radius) * z_fraction

            distance_map = _distance_to_membrane(
                R,
                abs_z,
                local_pore_radius,
                self.membrane_half_thickness,
            )
            conductivity_grid = _conductivity_from_distance(
                distance_map,
                self.bulk_conductivity,
                self.membrane_conductivity,
            )
            
            # Extract edges
            x_range = np.unique(self.X[:, 0, 0])
            y_range = np.unique(self.Y[0, :, 0])
            z_range = np.unique(self.Z[0, 0, :])
            
            self._interpolator = RegularGridInterpolator(
                (x_range, y_range, z_range), conductivity_grid,
                bounds_error=False, fill_value=self.bulk_conductivity
            )
        return self._interpolator

class BinFilePore(BasePore):
    def __init__(self, bin_file_path, base_sigma, mask_radius=-1, data_units="distance"):
        logger.info(f"Loading pore structure from binary file: {bin_file_path}")
        val3d, self.dimensions, self.grid_shape, metadata = readbinGrid(
            bin_file_path, mask_radius, return_metadata=True
        )
        self.grid_spacing = metadata["spacing"][0] if metadata else None
        self.base_sigma = base_sigma
        self.data_units = (data_units or "distance").lower()
        if self.data_units == "distance":
            calcSig = base_sigma * condfrac(val3d)
            logger.info("Binary file interpreted as distance map (condfrac × bulk conductivity).")
        elif self.data_units == "conductivity":
            calcSig = val3d.astype(np.float32)
            logger.info("Binary file interpreted as conductivity map (values used directly).")
        else:
            raise ValueError(f"Unsupported data_units '{self.data_units}' for BinFilePore.")
        grid_axes = (
            np.linspace(-self.dimensions[0]/2., self.dimensions[0]/2., num=self.grid_shape[0]),
            np.linspace(-self.dimensions[1]/2., self.dimensions[1]/2., num=self.grid_shape[1]),
            np.linspace(-self.dimensions[2]/2., self.dimensions[2]/2., num=self.grid_shape[2]),
        )
        self.interp = RegularGridInterpolator(
            grid_axes,
            calcSig,
            bounds_error=False,
            fill_value=base_sigma,
        )
        if self.data_units == "distance":
            self.distance_interp = RegularGridInterpolator(
                grid_axes,
                val3d,
                bounds_error=False,
                fill_value=0.0,
            )
        else:
            # Approximate distance map for radius checks using condfrac inversion.
            minr = 1.3
            maxr = 4.1
            slope = 1.0 / (maxr - minr)
            int_val = -minr * slope
            fraction = np.clip(calcSig / base_sigma, 0.0, 1.0)
            approx_distance = (fraction - int_val) / slope
            approx_distance = np.maximum(approx_distance, 0.0)
            self.distance_interp = RegularGridInterpolator(
                grid_axes,
                approx_distance,
                bounds_error=False,
                fill_value=0.0,
            )
            logger.info("Derived approximate distance map from conductivity values for radius checks.")
        logger.info(f"Binary file pore loaded successfully")
        logger.info(f"Grid dimensions: {self.grid_shape[0]} x {self.grid_shape[1]} x {self.grid_shape[2]}")
        logger.info(f"Effective size: {self.dimensions[0]:.1f} x {self.dimensions[1]:.1f} x {self.dimensions[2]:.1f} Å")
        logger.info(f"Conductivity range: {np.min(calcSig):.6f} to {np.max(calcSig):.6f} S/m")
    
    def get_conductivity_interpolator(self):
        return self.interp
    
    def get_dimensions(self):
        return self.dimensions
    
    def get_grid_shape(self):
        return self.grid_shape
    
    def get_grid_spacing(self):
        return self.grid_spacing

    def get_distance_interpolator(self):
        return self.distance_interp

class BiologicalPore(BasePore):
    def __init__(self, X, Y, Z, pore_pdb, membrane_half_thickness, 
                 bulk_conductivity=10.5, cutoff=5.0, 
                 use_vdw_radii=True, default_radius=1.5,
                 membrane_conductivity=0.0001, membrane_z_offset=0.0,
                 use_charges=False, debye_length=2.15, bjerrum_length=7.15,
                 resolution=1.0, cleanup_temp_files=True, box_dimensions=None,
                 temp_file_prefix="biological_pore",
                 use_direct_distance_calculation=False,
                 use_pdb2pqr=False, force_field='CHARMM', ph=7.0):
        try:
            import MDAnalysis as mda
        except Exception as exc:
            raise ImportError(
                "MDAnalysis is required to construct a BiologicalPore."
            ) from exc
        import subprocess
        import os
        from pathlib import Path
        import time
        
        self.X = X
        self.Y = Y
        self.Z = Z
        self.bulk_conductivity = bulk_conductivity
        self.membrane_conductivity = membrane_conductivity
        self._interpolator = None
        self.phi_interp = lambda pts: np.zeros(len(pts))
        self.dimensions = None
        self.grid_shape = None
        
        logger.info(f"Creating biological pore from {pore_pdb}")
        logger.info(f"Method: {'Direct distance calculation' if use_direct_distance_calculation else 'Subprocess PDB→XYZ→BIN'}")
        
        prepared_pore: Optional[PreparedStructure] = None
        pore_file = pore_pdb
        if use_pdb2pqr:
            logger.info(
                "Preparing pore structure with external pdb2pqr pipeline "
                "and custom radius overrides."
            )
            try:
                prepared_pore = prepare_structure(
                    pore_pdb,
                    ph=ph,
                    default_radius=default_radius,
                    use_external_pdb2pqr=True,
                    pdb2pqr_force_field=force_field,
                )
                pore_file = prepared_pore.pqr_file
            except ImportError as exc:
                logger.error(
                    "Structure preparation requires pdbfixer and openmm: %s", exc
                )
                raise
            except Exception as exc:
                logger.error(f"Failed to prepare pore structure: {exc}")
                raise
       
        try:
            pore_universe = mda.Universe(str(pore_file))
            pore_atoms = pore_universe.atoms
            logger.info(f"Loaded pore with {len(pore_atoms)} atoms")
        except Exception as e:
            logger.error(f"Failed to load pore file {pore_file}: {e}")
            raise
        
        if use_pdb2pqr:
            try:
                pore_radii = pore_atoms.radii.copy()
                if use_charges:
                    pore_charges = pore_atoms.charges.copy()
                    if not np.any(pore_charges):
                        logger.warning(
                            "Prepared pore PQR contains zero charges; charge-aware calculations "
                            "will assume neutral atoms."
                        )
                else:
                    pore_charges = np.zeros(len(pore_atoms))
                logger.info("Using radii from prepared PQR file")
            except Exception as exc:
                logger.error(f"Failed to read radii from prepared PQR: {exc}")
                raise
        elif use_vdw_radii:
            logger.info("Assigning van der Waals radii to pore atoms...")
            pore_radii = VanDerWaalsRadii.assign_radii_to_atoms(
                pore_atoms, 
                default_radius=default_radius,
                verbose=True
            )
            pore_charges = np.zeros(len(pore_atoms))
        else:
            try:
                pore_radii = pore_atoms.radii
                pore_charges = np.zeros(len(pore_atoms))
                logger.info("Using radii from pore PDB file")
            except:
                pore_radii = np.ones(len(pore_atoms)) * default_radius
                pore_charges = np.zeros(len(pore_atoms))
                logger.info(f"Using default radius {default_radius} Å for all pore atoms")
        
        pore_positions = pore_atoms.positions
        logger.info(f"Pore center of mass: {pore_atoms.center_of_mass()}")
        logger.info(f"Pore extent: X[{np.min(pore_positions[:,0]):.1f}, {np.max(pore_positions[:,0]):.1f}], "
                    f"Y[{np.min(pore_positions[:,1]):.1f}, {np.max(pore_positions[:,1]):.1f}], "
                    f"Z[{np.min(pore_positions[:,2]):.1f}, {np.max(pore_positions[:,2]):.1f}]")
        self.pore_positions = pore_positions
        self.pore_radii = pore_radii
        self.pore_tree = KDTree(pore_positions) if len(pore_positions) > 0 else None
        
        if membrane_half_thickness == 0.0:
            logger.info("Zero membrane thickness - creating interpolator using distance field")
            
            if box_dimensions is not None:
                x_min, x_max = box_dimensions['x']
                y_min, y_max = box_dimensions['y']
                z_min, z_max = box_dimensions['z']
            else:
                x_min, x_max = np.min(X), np.max(X)
                y_min, y_max = np.min(Y), np.max(Y)
                z_min, z_max = np.min(Z), np.max(Z)
            
            nx = int(np.ceil((x_max - x_min) / resolution))
            ny = int(np.ceil((y_max - y_min) / resolution))
            nz = int(np.ceil((z_max - z_min) / resolution))
            
            Lm = nx * resolution
            Wm = ny * resolution  
            Hm = nz * resolution
            
            self.dimensions = [Lm, Wm, Hm]
            self.grid_shape = [nx, ny, nz]
            
            if use_direct_distance_calculation:
                # Direct distance calculation approach
                logger.info("Using direct distance calculation from MDAnalysis...")
                
                # Create coordinate grids
                x_grid = np.linspace(-Lm/2., Lm/2., num=nx)
                y_grid = np.linspace(-Wm/2., Wm/2., num=ny)
                z_grid = np.linspace(-Hm/2., Hm/2., num=nz)
                
                logger.info("Creating coordinate grid...")
                X_val3d, Y_val3d, Z_val3d = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')
                grid_coords = np.column_stack([X_val3d.ravel(), Y_val3d.ravel(), Z_val3d.ravel()])
                
                logger.info(f"Computing distance field for {len(grid_coords):,} grid points using gen_dist.py-style algorithm...")

                # Match gen_dist.py parameters exactly
                SEARCH_LENGTH = 2  # Same as gen_dist.py
                max_search_radius = cutoff  # User's modification
                
                # Apply translation to match gen_dist.py workflow.
                # gen_dist.py translates atoms by subtracting the lower-bound
                # corner so they live in [0, Lm] x [0, Wm] x [0, Hm], which is
                # the coordinate space the spatial-hash indexing below assumes.
                # (Variable name predates this comment; kept for minimal diff.)
                center_coords = np.array([x_min, y_min, z_min])
                translated_pore_positions = pore_positions - center_coords
                
                # Set up spatial hashing like gen_dist.py
                cell_size = cutoff / SEARCH_LENGTH
                box_x, box_y, box_z = Lm, Wm, Hm
                cell_num_x = int(box_x * SEARCH_LENGTH / cutoff)
                cell_num_y = int(box_y * SEARCH_LENGTH / cutoff) 
                cell_num_z = int(box_z * SEARCH_LENGTH / cutoff)
                
                logger.info(f"Spatial hash parameters:")
                logger.info(f"  SEARCH_LENGTH: {SEARCH_LENGTH}")
                logger.info(f"  cell_size: {cell_size:.3f}")
                logger.info(f"  cells: {cell_num_x} x {cell_num_y} x {cell_num_z}")
                
                def get_cell_index(x, y, z):
                    """Get cell index for spatial hashing (matches gen_dist.py)"""
                    cell_x = min(int(x / cell_size), cell_num_x - 1) if cell_num_x > 0 else 0
                    cell_y = min(int(y / cell_size), cell_num_y - 1) if cell_num_y > 0 else 0
                    cell_z = min(int(z / cell_size), cell_num_z - 1) if cell_num_z > 0 else 0
                    return cell_x + cell_num_x * cell_y + cell_num_x * cell_num_y * cell_z
                
                # Build spatial hash (same as gen_dist.py)
                logger.info("Building spatial hash...")
                cell_atoms = {}
                for i, pos in enumerate(translated_pore_positions):
                    cell_idx = get_cell_index(pos[0], pos[1], pos[2])
                    if cell_idx not in cell_atoms:
                        cell_atoms[cell_idx] = []
                    cell_atoms[cell_idx].append(i)
                
                logger.info(f"Created {len(cell_atoms)} non-empty cells")
                
                # Process grid points using same logic as gen_dist.py
                val3d_flat = np.full(len(grid_coords), cutoff, dtype=np.float32)
                
                for grid_idx, grid_coord in enumerate(grid_coords):
                    # Translate grid coordinate to match gen_dist.py coordinate system
                    translated_coord = grid_coord  # Already centered around origin
                    
                    # Find which cell this grid point belongs to
                    cell_loc_x = min(int(translated_coord[0] / cell_size), cell_num_x - 1) if cell_num_x > 0 else 0
                    cell_loc_y = min(int(translated_coord[1] / cell_size), cell_num_y - 1) if cell_num_y > 0 else 0
                    cell_loc_z = min(int(translated_coord[2] / cell_size), cell_num_z - 1) if cell_num_z > 0 else 0
                    
                    # Search neighboring cells (same range as gen_dist.py)
                    start_x = max(0, cell_loc_x - SEARCH_LENGTH)
                    end_x = min(cell_num_x, cell_loc_x + SEARCH_LENGTH + 1)
                    start_y = max(0, cell_loc_y - SEARCH_LENGTH)
                    end_y = min(cell_num_y, cell_loc_y + SEARCH_LENGTH + 1)
                    start_z = max(0, cell_loc_z - SEARCH_LENGTH)
                    end_z = min(cell_num_z, cell_loc_z + SEARCH_LENGTH + 1)
                    
                    min_dist = cutoff
                    
                    # Search all neighboring cells
                    for z in range(start_z, end_z):
                        for y in range(start_y, end_y):
                            for x in range(start_x, end_x):
                                cell_hash = x + cell_num_x * y + cell_num_x * cell_num_y * z
                                
                                if cell_hash in cell_atoms:
                                    for atom_idx in cell_atoms[cell_hash]:
                                        atom_pos = translated_pore_positions[atom_idx]
                                        atom_radius = pore_radii[atom_idx]
                                        
                                        # Calculate distance (same as gen_dist.py)
                                        dx = atom_pos[0] - translated_coord[0]
                                        dy = atom_pos[1] - translated_coord[1]
                                        dz = atom_pos[2] - translated_coord[2]
                                        
                                        dist_to_center = np.sqrt(dx*dx + dy*dy + dz*dz)
                                        dist_to_surface = dist_to_center - atom_radius
                                        
                                        min_dist = min(min_dist, dist_to_surface)
                    
                    # Store result (ensure non-negative, same as gen_dist.py)
                    val3d_flat[grid_idx] = max(min_dist, 0.0)
                    
                    if grid_idx % 100000 == 0 and grid_idx > 0:
                        logger.debug(f"  Processed {grid_idx:,} / {len(grid_coords):,} points")
                
                # Reshape to 3D (same order as gen_dist.py)
                val3d = val3d_flat.reshape((nx, ny, nz))
                
                # Apply cutoff clamping like gen_dist.py
                val3d = np.minimum(val3d, cutoff).astype(np.float32)
                
                logger.info(f"Direct distance calculation complete!")
                logger.info(f"Distance field range: [{np.min(val3d):.3f}, {np.max(val3d):.3f}] Å")
                logger.info(f"Distance field shape: {val3d.shape}")
                logger.info(f"Using gen_dist.py-compatible algorithm")
                
            else:
                # Subprocess approach for val3d creation (PDB → XYZ → BIN)
                logger.info("Using subprocess PDB→XYZ→BIN approach...")
                
                # Find scripts (same as before)
                try:
                    from sem.scripts import PDB2XYZ_SCRIPT, GEN_DIST_SCRIPT
                    pdb2xyz_script = PDB2XYZ_SCRIPT
                    gen_dist_script = GEN_DIST_SCRIPT
                    logger.info(f"Using scripts from installed package")
                except ImportError:
                    script_locations = [
                        Path(__file__).parent / "scripts",
                        Path(__file__).parent,
                        Path.cwd(),
                        Path("scripts"),
                    ]
                    
                    pdb2xyz_script = None
                    gen_dist_script = None
                    
                    for script_dir in script_locations:
                        potential_pdb2xyz = script_dir / "pdb2xyz.py"
                        potential_gen_dist = script_dir / "gen_dist.py"
                        
                        if potential_pdb2xyz.exists() and potential_gen_dist.exists():
                            pdb2xyz_script = potential_pdb2xyz
                            gen_dist_script = potential_gen_dist
                            logger.info(f"Found scripts in: {script_dir}")
                            break
                    
                    if not pdb2xyz_script or not gen_dist_script:
                        raise FileNotFoundError(
                            "pdb2xyz.py and gen_dist.py not found. "
                            "Please ensure they are installed with the package or available in the current directory."
                        )
                
                # Create temporary files
                timestamp = int(time.time())
                process_id = os.getpid()
                xyz_file = Path(f"{temp_file_prefix}_{timestamp}_{process_id}.xyz")
                bin_file = Path(f"{temp_file_prefix}_{timestamp}_{process_id}.bin")
                
                logger.info(f"Using local temporary files:")
                logger.info(f"   XYZ file: {xyz_file}")
                logger.info(f"   BIN file: {bin_file}")
                
                try:
                    # Step 1: Convert to XYZ using alternative process if use_pdb2pqr
                    logger.info("Step 1: Converting to XYZ format...")
                    if use_pdb2pqr and prepared_pore is not None:
                        # Alternative: Use in-memory PQR to XYZ conversion instead of script
                        logger.info("Using alternative PQR to XYZ conversion")
                        universe = mda.Universe(str(prepared_pore.pqr_file))
                        selection = 'not resname HOH WAT TIP3 SOL NA CL K MG CA ZN'
                        atoms = universe.select_atoms(selection)
                        
                        positions = atoms.positions
                        radii = atoms.radii  # Radii from B-factor
                        
                        with open(xyz_file, 'w') as f:
                            for i in range(len(atoms)):
                                x, y, z = positions[i]
                                r = radii[i]
                                line = f"{x:.3f} {y:.3f} {z:.3f} {r:.3f}"
                                f.write(line + '\n')
                        
                        logger.info(f"Converted PQR to XYZ: {len(atoms)} atoms")
                    else:
                        # Original: Use pdb2xyz script on PDB
                        pdb2xyz_cmd = [sys.executable, str(pdb2xyz_script), str(pore_pdb), str(xyz_file)]
                        result = subprocess.run(pdb2xyz_cmd, capture_output=True, text=True, timeout=300)
                        
                        if result.returncode != 0:
                            logger.error(f"pdb2xyz.py failed with return code {result.returncode}")
                            logger.error(f"STDERR: {result.stderr}")
                            raise RuntimeError("PDB to XYZ conversion failed")
                    
                    # Step 2: XYZ → BIN
                    logger.info("Step 2: Converting XYZ to binary distance field...")
                    gen_dist_cmd = [
                        sys.executable, str(gen_dist_script), str(xyz_file),
                        str(x_max), str(y_max), str(z_max),
                        str(x_min), str(y_min), str(z_min),
                        str(resolution), str(cutoff), str(bin_file)
                    ]
                    result = subprocess.run(gen_dist_cmd, capture_output=True, text=True, timeout=600)
                    
                    if result.returncode != 0:
                        logger.error(f"gen_dist.py failed with return code {result.returncode}")
                        logger.error(f"STDERR: {result.stderr}")
                        raise RuntimeError("XYZ to binary distance field conversion failed")
                    
                    # Step 3: Read BIN → val3d
                    logger.info("Step 3: Reading binary file to get val3d...")
                    val3d, [Lm, Wm, Hm], [nx, ny, nz] = readbinGrid(str(bin_file), mask_radius=-1)
                    
                    logger.info(f"Subprocess approach complete!")
                    logger.info(f"Distance field range: [{np.min(val3d):.3f}, {np.max(val3d):.3f}] Å")
                    
                finally:
                    # Cleanup temporary files
                    if cleanup_temp_files:
                        try:
                            if xyz_file.exists():
                                xyz_file.unlink()
                            if bin_file.exists():
                                bin_file.unlink()
                            if use_pdb2pqr and prepared_pore is not None:
                                prepared_pore.cleanup()
                        except Exception as cleanup_error:
                            logger.warning(f"Failed to cleanup temporary files: {cleanup_error}")
                
            calcSig = bulk_conductivity * condfrac(val3d)
            self._interpolator = RegularGridInterpolator(
                (np.linspace(-Lm/2., Lm/2., num=nx),
                 np.linspace(-Wm/2., Wm/2., num=ny),
                 np.linspace(-Hm/2., Hm/2., num=nz)),
                calcSig,
                bounds_error=False,
                fill_value=bulk_conductivity
            )
            
        else:
            # Case 2: With membrane - use cylindrical-style approach
            logger.info(f"Creating biological pore with membrane (thickness={membrane_half_thickness*2:.1f}Å)")
            logger.info("Using cylindrical-style mask creation approach")
            
            # Calculate R coordinate (radial distance from Z axis) - same as cylindrical
            R = np.sqrt(X**2 + Y**2)
            
            # Initialize with bulk conductivity
            base_conductivity = np.ones(X.shape) * bulk_conductivity
            
            # Calculate displaced Z coordinates for membrane positioning
            Z_displaced = Z - membrane_z_offset
            
            # Calculate local pore radius for each Z slice
            logger.info("Calculating local pore radius based on biological structure...")
            
            z_coords = Z[0, 0, :]  # Get Z coordinates from simulation grid
            local_pore_radius = np.zeros_like(z_coords)
            
            # For each Z slice, find the maximum extent of pore atoms
            for i, z_pos in enumerate(z_coords):
                # Check if this z position is within the displaced membrane
                if np.abs(z_pos - membrane_z_offset) < membrane_half_thickness:
                    # Find pore atoms near this Z position (within ±2 Å slice)
                    z_slice_mask = np.abs(pore_positions[:, 2] - z_pos) <= 2.0
                    
                    if np.any(z_slice_mask):
                        # Get atoms in this slice
                        slice_positions = pore_positions[z_slice_mask]
                        slice_radii = pore_radii[z_slice_mask]
                        
                        # Calculate radial distance from Z-axis for each atom
                        atom_radial_distances = np.sqrt(slice_positions[:, 0]**2 + slice_positions[:, 1]**2)
                        
                        # Maximum extent = radial distance + atom radius
                        max_extents = atom_radial_distances + slice_radii
                        
                        # Use the minimum extent as local pore radius
                        local_pore_radius[i] = np.min(max_extents) if len(max_extents) > 0 else 0
                    else:
                        # No atoms in this slice
                        local_pore_radius[i] = 0
                else:
                    # Outside membrane region - mark for later interpolation
                    local_pore_radius[i] = np.nan

            # Propagate nearest valid radius outside the membrane region so the pore column
            # remains well-defined for transition calculations.
            radius_indices = np.arange(len(local_pore_radius))
            valid_mask = np.isfinite(local_pore_radius)
            if np.any(valid_mask):
                local_pore_radius = np.interp(
                    radius_indices,
                    radius_indices[valid_mask],
                    local_pore_radius[valid_mask]
                )
            else:
                box_max_radius = min(np.max(X) - np.min(X), np.max(Y) - np.min(Y)) / 2
                local_pore_radius[:] = box_max_radius
            
            # Create 3D local_pore_radius array matching grid shape (same as cylindrical approach)
            local_pore_radius_3d = np.zeros_like(Z)
            for i, z_pos in enumerate(z_coords):
                local_pore_radius_3d[:, :, i] = local_pore_radius[i]
            
            # Create membrane mask using displaced coordinates (EXACTLY like cylindrical approach)
            # Membrane exists where: within displaced membrane thickness AND outside local pore radius
            membrane_mask = (np.abs(Z_displaced) < membrane_half_thickness) & (R > local_pore_radius_3d)
            
            # Apply conductivity model for fine-grained pore structure within the allowed pore region
            if len(pore_atoms) > 0:
                logger.info("Applying fine-grained biological pore structure...")
                
                # Create conductivity model for pore
                conductivity_model = SimpleConductivityModel(
                    bulk_conductivity=bulk_conductivity,
                )
                
                # Flatten grid coordinates for processing
                grid_coords = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
                
                # Initialize with bulk conductivity
                conductivity_map = np.full(len(grid_coords), bulk_conductivity)
                
                # Only modify conductivity within the allowed pore region (not in membrane)
                pore_region_flat = ~membrane_mask.ravel()
                
                if np.any(pore_region_flat):
                    pore_tree = KDTree(pore_positions)
                    distances, indices = pore_tree.query(
                        grid_coords[pore_region_flat], 
                        distance_upper_bound=cutoff + np.max(pore_radii)
                    )
                    
                    # Adjust distances for atom radii and calculate conductivity
                    valid_mask = ~np.isinf(distances)
                    if np.any(valid_mask):
                        distances[valid_mask] -= pore_radii[indices[valid_mask]]
                        distances[valid_mask] = np.maximum(distances[valid_mask], 0)
                        
                        # Apply conductivity model
                        modified_cond = conductivity_model(distances[valid_mask])
                        
                        # Update conductivity only for valid pore region points
                        pore_indices = np.where(pore_region_flat)[0][valid_mask]
                        conductivity_map[pore_indices] = np.minimum(
                            conductivity_map[pore_indices],
                            modified_cond
                        )
                
                # Reshape back to grid shape
                conductivity_map = conductivity_map.reshape(X.shape)
            else:
                conductivity_map = base_conductivity.copy()
            
            # Apply smooth membrane → bulk transition outside the biological pore interior
            distance_map = _distance_to_membrane(
                R,
                np.abs(Z_displaced),
                local_pore_radius_3d,
                membrane_half_thickness,
            )
            transition_conductivity = _conductivity_from_distance(
                distance_map,
                bulk_conductivity,
                membrane_conductivity,
            )
            pore_region_mask = (np.abs(Z_displaced) < membrane_half_thickness) & (R <= local_pore_radius_3d)
            apply_transition_mask = ~pore_region_mask
            conductivity_map[apply_transition_mask] = np.minimum(
                conductivity_map[apply_transition_mask],
                transition_conductivity[apply_transition_mask],
            )
            
            # Log statistics
            n_total = X.size
            n_membrane = np.sum(membrane_mask)
            n_pore = np.sum(~membrane_mask & (np.abs(Z_displaced) < membrane_half_thickness))
            n_bulk = n_total - n_membrane - n_pore
            
            logger.info(f"Biological pore conductivity map statistics:")
            logger.info(f"  Total grid points: {n_total:,}")
            logger.info(f"  Membrane: {n_membrane:,} points ({100*n_membrane/n_total:.1f}%)")
            logger.info(f"  Pore: {n_pore:,} points ({100*n_pore/n_total:.1f}%)")
            logger.info(f"  Bulk: {n_bulk:,} points ({100*n_bulk/n_total:.1f}%)")
            logger.info(f"  Max local pore radius: {np.max(local_pore_radius):.1f} Å")
            logger.info(f"  Min local pore radius: {np.min(local_pore_radius[local_pore_radius > 0]):.1f} Å")
            logger.info(f"  Membrane z-offset: {membrane_z_offset:.1f} Å")
            logger.info(f"  Membrane placed outside biological pore structure")
            
            self.conductivity_grid = conductivity_map  # Assume from logic
            # self.phi_interp = phi_interp
            
            # Set dimensions from box
            if box_dimensions:
                self.dimensions = [
                    box_dimensions['x'][1] - box_dimensions['x'][0],
                    box_dimensions['y'][1] - box_dimensions['y'][0],
                    box_dimensions['z'][1] - box_dimensions['z'][0]
                ]
                self.grid_shape = self.X.shape

        if cleanup_temp_files and use_pdb2pqr and prepared_pore is not None:
            prepared_pore.cleanup()
    
    def get_conductivity_interpolator(self):
        if self._interpolator is None:
            # Create from conductivity_grid
            x_range = np.unique(self.X[:, 0, 0])
            y_range = np.unique(self.Y[0, :, 0])
            z_range = np.unique(self.Z[0, 0, :])
            self._interpolator = RegularGridInterpolator(
                (x_range, y_range, z_range), self.conductivity_grid,
                bounds_error=False, fill_value=self.bulk_conductivity
            )
        return self._interpolator
    
    def get_phi_interpolator(self):
        return self.phi_interp
    
    def get_dimensions(self):
        return self.dimensions
    
    def get_grid_shape(self):
        return self.grid_shape

class PoreGeometry:
    @classmethod
    def create_pore(cls, pore_type, X=None, Y=None, Z=None, **kwargs):
        pore_type = pore_type.lower()
        if pore_type == "cylindrical":
            return CylindricalPore(X, Y, Z, **kwargs)
        elif pore_type == "double_cone":
            return DoubleConePore(X, Y, Z, **kwargs)
        elif pore_type == "bin_file":
            return BinFilePore(**kwargs)
        elif pore_type == "biological":
            return BiologicalPore(X, Y, Z, **kwargs)
        raise ValueError(f"Unknown pore_type: {pore_type}")
