#!/usr/bin/env python3
"""
Python implementation of gen_dist for generating distance fields from XYZ files.
Equivalent to the C version for cross-platform compatibility.

Usage: python gen_dist.py <xyz_file> <MinX> <MinY> <MinZ> <MaxX> <MaxY> <MaxZ> <Resolution> <cutoff> <OutputFile>
"""

import sys
import numpy as np
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from numba import jit, prange, types
from numba.typed import List

class Atom:
    """Atom structure - kept as simple Python class for interface compatibility"""
    def __init__(self, x, y, z, radius):
        self.x = float(x)
        self.y = float(y) 
        self.z = float(z)
        self.radius = float(radius)
        self.hVal = 0

# Numba-compiled functions for critical path optimization
@jit(nopython=True, fastmath=False, parallel=False)
def compute_cell_index(x, y, z, cell_size, cell_num_x, cell_num_y, cell_num_z):
    """Compute cell index - exact same logic as original, but compiled"""
    cell_x = min(int(x / cell_size), cell_num_x - 1)
    cell_y = min(int(y / cell_size), cell_num_y - 1) 
    cell_z = min(int(z / cell_size), cell_num_z - 1)
    
    cell_x = max(0, cell_x)
    cell_y = max(0, cell_y)
    cell_z = max(0, cell_z)
    
    return cell_x + cell_num_x * cell_y + cell_num_x * cell_num_y * cell_z

@jit(nopython=True, fastmath=False)
def compute_distance_slice_numba(
    grid_point_num_x, grid_point_num_y, resolution, cutoff,
    cell_size, cell_num_x, cell_num_y, cell_num_z, SEARCH_LENGTH,
    z_start, z_end,
    atom_positions, atom_radii, atom_hash_values, cell_starts, cell_ends
):
    """
    Numba-compiled slice computation that maintains EXACT same logic as original.
    No approximations or algorithmic changes - just compiled for speed.
    """
    slice_size = z_end - z_start
    result = np.full((slice_size, grid_point_num_y, grid_point_num_x), cutoff, dtype=np.float32)
    
    for grid_z_idx in range(z_start, z_end):
        grid_coord_z = grid_z_idx * resolution
        
        # Exact same cell calculation logic as original
        cell_loc_z = min(int(grid_coord_z / cell_size), cell_num_z - 1)
        cell_loc_z = max(0, cell_loc_z)
        start_z = max(0, cell_loc_z - SEARCH_LENGTH)
        end_z = min(cell_num_z, cell_loc_z + SEARCH_LENGTH + 1)
        
        for grid_y_idx in range(grid_point_num_y):
            grid_coord_y = grid_y_idx * resolution
            cell_loc_y = min(int(grid_coord_y / cell_size), cell_num_y - 1)
            cell_loc_y = max(0, cell_loc_y)
            start_y = max(0, cell_loc_y - SEARCH_LENGTH)
            end_y = min(cell_num_y, cell_loc_y + SEARCH_LENGTH + 1)
            
            for grid_x_idx in range(grid_point_num_x):
                grid_coord_x = grid_x_idx * resolution
                cell_loc_x = min(int(grid_coord_x / cell_size), cell_num_x - 1)
                cell_loc_x = max(0, cell_loc_x)
                start_x = max(0, cell_loc_x - SEARCH_LENGTH)
                end_x = min(cell_num_x, cell_loc_x + SEARCH_LENGTH + 1)
                
                min_dist = cutoff
                
                # Triple nested loop - exact same as original
                for z in range(start_z, end_z):
                    for y in range(start_y, end_y):
                        for x in range(start_x, end_x):
                            # Same hash calculation
                            cell_hash = x + cell_num_x * y + cell_num_x * cell_num_y * z
                            
                            # Check if this cell has atoms
                            if cell_hash < len(cell_starts) and cell_starts[cell_hash] >= 0:
                                start_atom = cell_starts[cell_hash]
                                end_atom = cell_ends[cell_hash]
                                
                                # Check each atom in cell - exact same calculation
                                for atom_idx in range(start_atom, end_atom):
                                    atom_x = atom_positions[atom_idx, 0]
                                    atom_y = atom_positions[atom_idx, 1]
                                    atom_z = atom_positions[atom_idx, 2]
                                    atom_radius = atom_radii[atom_idx]
                                    
                                    # Exact same distance calculation
                                    dx = atom_x - grid_coord_x
                                    dy = atom_y - grid_coord_y
                                    dz = atom_z - grid_coord_z
                                    
                                    dist_to_center = np.sqrt(dx*dx + dy*dy + dz*dz)
                                    dist_to_surface = dist_to_center - atom_radius
                                    
                                    if dist_to_surface < min_dist:
                                        min_dist = dist_to_surface
                
                # Store result - exact same logic
                result[grid_z_idx - z_start, grid_y_idx, grid_x_idx] = max(min_dist, 0.0)
    
    return result

class PrecisionOptimizedGenerator:
    """
    Optimized generator that maintains exact mathematical precision.
    Only speeds up computation without changing any algorithms.
    """
    
    def __init__(self, resolution, cutoff, x_size, y_size, z_size):
        self.resolution = resolution
        self.cutoff = cutoff
        self.x_size = x_size
        self.y_size = y_size
        self.z_size = z_size
        self.SEARCH_LENGTH = 2  # Exact same as original
        
        # Exact same calculations as original
        self.cell_size = cutoff / self.SEARCH_LENGTH
        self.cell_num_x = int(x_size * self.SEARCH_LENGTH / cutoff)
        self.cell_num_y = int(y_size * self.SEARCH_LENGTH / cutoff)
        self.cell_num_z = int(z_size * self.SEARCH_LENGTH / cutoff)
        
        self.grid_point_num_x = int(x_size / resolution) + 1
        self.grid_point_num_y = int(y_size / resolution) + 1
        self.grid_point_num_z = int(z_size / resolution) + 1
        
        # Minimal output like original
        pass
    
    def _convert_to_numba_format(self, atoms):
        """Convert atoms to numba-friendly arrays while preserving exact data"""
        n_atoms = len(atoms)
        
        # Pre-allocate arrays
        atom_positions = np.zeros((n_atoms, 3), dtype=np.float64)  # Use float64 for precision
        atom_radii = np.zeros(n_atoms, dtype=np.float64)
        atom_hash_values = np.zeros(n_atoms, dtype=np.int32)
        
        # Fill arrays with exact same data
        for i, atom in enumerate(atoms):
            atom_positions[i, 0] = atom.x
            atom_positions[i, 1] = atom.y
            atom_positions[i, 2] = atom.z
            atom_radii[i] = atom.radius
            atom_hash_values[i] = atom.hVal
        
        return atom_positions, atom_radii, atom_hash_values
    
    def _setup_spatial_hash_optimized(self, atoms):
        """Optimized spatial hash setup maintaining exact same logic"""
        # Minimal output - no debug prints
        
        # Assign hash values - exact same logic as original
        for atom in atoms:
            atom.hVal = compute_cell_index(
                atom.x, atom.y, atom.z, self.cell_size, 
                self.cell_num_x, self.cell_num_y, self.cell_num_z
            )
        
        # Sort by hash value - exact same as original
        atoms.sort(key=lambda a: a.hVal)
        
        # Convert to numba format after sorting
        atom_positions, atom_radii, atom_hash_values = self._convert_to_numba_format(atoms)
        
        # Create cell lookup tables
        max_hash = self.cell_num_x * self.cell_num_y * self.cell_num_z
        cell_starts = np.full(max_hash, -1, dtype=np.int32)
        cell_ends = np.full(max_hash, -1, dtype=np.int32)
        
        if len(atoms) > 0:
            current_hash = atom_hash_values[0]
            start_idx = 0
            
            for i in range(len(atoms)):
                if atom_hash_values[i] != current_hash:
                    cell_starts[current_hash] = start_idx
                    cell_ends[current_hash] = i
                    current_hash = atom_hash_values[i]
                    start_idx = i
            
            # Handle last group
            cell_starts[current_hash] = start_idx
            cell_ends[current_hash] = len(atoms)
        
        return atom_positions, atom_radii, cell_starts, cell_ends
    
    def _compute_distance_slice_optimized(self, atom_data, z_start, z_end):
        """Optimized slice computation maintaining exact precision"""
        atom_positions, atom_radii, cell_starts, cell_ends = atom_data
        
        return compute_distance_slice_numba(
            self.grid_point_num_x, self.grid_point_num_y, self.resolution, self.cutoff,
            self.cell_size, self.cell_num_x, self.cell_num_y, self.cell_num_z, self.SEARCH_LENGTH,
            z_start, z_end,
            atom_positions, atom_radii, np.arange(len(atom_radii)), cell_starts, cell_ends
        )
    
    def generate_distance_field(self, atoms, num_threads=None):
        """Generate distance field with speed optimization but exact precision"""
        if num_threads is None:
            num_threads = min(12, multiprocessing.cpu_count())  # Use more threads
        
        # Setup spatial hash with optimization
        atom_data = self._setup_spatial_hash_optimized(atoms)
        
        # Use more aggressive chunking for better parallelization
        chunk_size = max(2, self.grid_point_num_z // (num_threads * 3))
        chunks = []
        for i in range(0, self.grid_point_num_z, chunk_size):
            chunks.append((i, min(i + chunk_size, self.grid_point_num_z)))
        
        start_time = time.time()
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            for z_start, z_end in chunks:
                future = executor.submit(self._compute_distance_slice_optimized, atom_data, z_start, z_end)
                futures.append((future, z_start, z_end))
            
            # Pre-allocate result array
            distance_field = np.zeros((self.grid_point_num_z, self.grid_point_num_y, self.grid_point_num_x), 
                                    dtype=np.float32)
            
            # Collect results efficiently
            for future, z_start, z_end in futures:
                result = future.result()
                distance_field[z_start:z_end] = result
        
        return distance_field

# Keep all the existing I/O functions exactly the same
def load_xyz_atoms(filename, x_lower, y_lower, z_lower, x_upper, y_upper, z_upper):
    """Load atoms - exact same logic as original"""
    atoms = []
    
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    parts = line.split()
                    if len(parts) >= 4:
                        x, y, z, r = map(float, parts[:4])
                        
                        if x > x_upper or y > y_upper or z > z_upper:
                            continue
                        
                        translated_x = x - x_lower
                        translated_y = y - y_lower  
                        translated_z = z - z_lower
                        
                        if translated_x > 0 and translated_y > 0 and translated_z > 0:
                            atom = Atom(x, y, z, r)
                            atoms.append(atom)
                            
                except ValueError as e:
                    continue
    
    except FileNotFoundError:
        print(f"ERROR: Cannot open file {filename}")
        sys.exit(1)
    
    return atoms

def apply_translation(atoms, x_lower, y_lower, z_lower):
    """Apply translation - exact same logic as original"""
    translated_atoms = []
    
    for atom in atoms:
        new_atom = Atom(
            atom.x - x_lower,
            atom.y - y_lower,
            atom.z - z_lower,
            atom.radius
        )
        translated_atoms.append(new_atom)
    
    return translated_atoms

def write_binary_file(distance_field, origin, resolution, filename):
    """Write binary file - exact same format as original"""
    try:
        with open(filename, 'wb') as f:
            grid_shape = distance_field.shape
            x_count = float(grid_shape[2])
            y_count = float(grid_shape[1])  
            z_count = float(grid_shape[0])
            
            f.write(np.array([x_count], dtype=np.float32).tobytes())
            f.write(np.array([y_count], dtype=np.float32).tobytes())
            f.write(np.array([z_count], dtype=np.float32).tobytes())
            f.write(np.array(origin, dtype=np.float32).tobytes())
            f.write(np.array([resolution], dtype=np.float32).tobytes())
            
            for z in range(grid_shape[0]):
                slice_data = distance_field[z, :, :].astype(np.float32)
                f.write(slice_data.tobytes())
        
    except Exception as e:
        print(f"ERROR: Failed to write {filename}: {e}")
        sys.exit(1)

def generate_binary_distance_field(
    xyz_file,
    x_lower,
    y_lower,
    z_lower,
    x_upper,
    y_upper,
    z_upper,
    resolution,
    cutoff,
    output_file,
    *,
    num_threads=None,
):
    """
    Programmatic helper that replicates the CLI workflow.
    Returns metadata about the generated grid so other scripts can inspect it.
    """
    xyz_file = Path(xyz_file)
    output_file = Path(output_file)

    box_x = x_upper - x_lower
    box_y = y_upper - y_lower
    box_z = z_upper - z_lower

    atoms = load_xyz_atoms(xyz_file, x_lower, y_lower, z_lower, x_upper, y_upper, z_upper)

    if not atoms:
        raise ValueError(f"No atoms loaded from {xyz_file}")
    
    translated_atoms = apply_translation(atoms, x_lower, y_lower, z_lower)
    
    generator = PrecisionOptimizedGenerator(resolution, cutoff + 2, box_x, box_y, box_z)
    distance_field = generator.generate_distance_field(translated_atoms, num_threads=num_threads)
    
    distance_field = np.minimum(distance_field, cutoff)
    
    origin = [x_lower, y_lower, z_lower]
    write_binary_file(distance_field, origin, resolution, output_file)

    metadata = {
        "output_file": str(output_file),
        "resolution": float(resolution),
        "origin": origin,
        "grid_shape": (
            int(distance_field.shape[2]),
            int(distance_field.shape[1]),
            int(distance_field.shape[0])
        ),
        "bounds": {
            "min": [x_lower, y_lower, z_lower],
            "max": [x_upper, y_upper, z_upper]
        }
    }
    return metadata

def main():
    """Main function - minimal output like original"""
    if len(sys.argv) != 11:
        print("Usage: python precision_optimized_gen_dist.py <xyz_file> <MinX> <MinY> <MinZ> <MaxX> <MaxY> <MaxZ> <Resolution> <cutoff> <OutputFile>")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    x_upper = float(sys.argv[2])
    y_upper = float(sys.argv[3])
    z_upper = float(sys.argv[4])
    x_lower = float(sys.argv[5])
    y_lower = float(sys.argv[6])
    z_lower = float(sys.argv[7])
    resolution = float(sys.argv[8])
    cutoff = float(sys.argv[9])
    output_file = sys.argv[10]
    
    try:
        generate_binary_distance_field(
            xyz_file,
            x_lower,
            y_lower,
            z_lower,
            x_upper,
            y_upper,
            z_upper,
            resolution,
            cutoff,
            output_file,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)

if __name__ == "__main__":
    main()
