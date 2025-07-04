#!/usr/bin/env python
"""
Modified SEM calculation for vertical movement of analyte through nanopore.
This version uses van der Waals radii for more accurate atomic volume calculations
and supports different pore geometries (cylindrical and double cone).
All units are in Angstroms throughout.

Enhanced with JSON configuration support and command-line interface.
"""

import os
import sys
import time
import json
import argparse
import logging
import numpy as np
from pathlib import Path
from fenics import *
import MDAnalysis as mda
from scipy.spatial import KDTree
from scipy.interpolate import RegularGridInterpolator

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s %(name)s: %(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

class VanDerWaalsRadii:
    """
    Van der Waals radii database for common elements in Angstroms.
    Based on Bondi (1964) and more recent literature values.
    """
    
    # Standard van der Waals radii in Angstroms
    VDW_RADII = {
        'H': 1.20,   # Hydrogen
        'HE': 1.40,  # Helium
        'LI': 1.82,  # Lithium
        'BE': 1.53,  # Beryllium
        'B': 1.92,   # Boron
        'C': 1.70,   # Carbon
        'N': 1.55,   # Nitrogen
        'O': 1.52,   # Oxygen
        'F': 1.47,   # Fluorine
        'NE': 1.54,  # Neon
        'NA': 2.27,  # Sodium
        'MG': 1.73,  # Magnesium
        'AL': 1.84,  # Aluminum
        'SI': 2.10,  # Silicon
        'P': 1.80,   # Phosphorus
        'S': 1.80,   # Sulfur
        'CL': 1.75,  # Chlorine
        'AR': 1.88,  # Argon
        'K': 2.75,   # Potassium
        'CA': 2.31,  # Calcium
        'SC': 2.11,  # Scandium
        'TI': 1.87,  # Titanium
        'V': 1.79,   # Vanadium
        'CR': 1.89,  # Chromium
        'MN': 1.97,  # Manganese
        'FE': 1.94,  # Iron
        'CO': 1.92,  # Cobalt
        'NI': 1.84,  # Nickel
        'CU': 1.32,  # Copper
        'ZN': 1.22,  # Zinc
        'GA': 1.87,  # Gallium
        'GE': 2.11,  # Germanium
        'AS': 1.85,  # Arsenic
        'SE': 1.90,  # Selenium
        'BR': 1.85,  # Bromine
        'KR': 2.02,  # Krypton
        'RB': 3.03,  # Rubidium
        'SR': 2.49,  # Strontium
        'Y': 2.32,   # Yttrium
        'ZR': 2.23,  # Zirconium
        'NB': 2.18,  # Niobium
        'MO': 2.17,  # Molybdenum
        'TC': 2.16,  # Technetium
        'RU': 2.13,  # Ruthenium
        'RH': 2.10,  # Rhodium
        'PD': 2.10,  # Palladium
        'AG': 1.72,  # Silver
        'CD': 1.58,  # Cadmium
        'IN': 1.93,  # Indium
        'SN': 2.17,  # Tin
        'SB': 2.06,  # Antimony
        'TE': 2.06,  # Tellurium
        'I': 1.98,   # Iodine
        'XE': 2.16,  # Xenon
        'CS': 3.43,  # Cesium
        'BA': 2.68,  # Barium
        'LA': 2.43,  # Lanthanum
        'CE': 2.42,  # Cerium
        'PR': 2.40,  # Praseodymium
        'ND': 2.39,  # Neodymium
        'PM': 2.38,  # Promethium
        'SM': 2.36,  # Samarium
        'EU': 2.35,  # Europium
        'GD': 2.34,  # Gadolinium
        'TB': 2.33,  # Terbium
        'DY': 2.31,  # Dysprosium
        'HO': 2.30,  # Holmium
        'ER': 2.29,  # Erbium
        'TM': 2.27,  # Thulium
        'YB': 2.26,  # Ytterbium
        'LU': 2.24,  # Lutetium
        'HF': 2.23,  # Hafnium
        'TA': 2.22,  # Tantalum
        'W': 2.18,   # Tungsten
        'RE': 2.16,  # Rhenium
        'OS': 2.16,  # Osmium
        'IR': 2.13,  # Iridium
        'PT': 2.13,  # Platinum
        'AU': 1.66,  # Gold
        'HG': 1.55,  # Mercury
        'TL': 1.96,  # Thallium
        'PB': 2.02,  # Lead
        'BI': 2.07,  # Bismuth
        'PO': 1.97,  # Polonium
        'AT': 2.02,  # Astatine
        'RN': 2.20,  # Radon
        'FR': 3.48,  # Francium
        'RA': 2.83,  # Radium
        'AC': 2.47,  # Actinium
        'TH': 2.45,  # Thorium
        'PA': 2.43,  # Protactinium
        'U': 2.41,   # Uranium
        'NP': 2.39,  # Neptunium
        'PU': 2.43,  # Plutonium
        'AM': 2.44,  # Americium
        'CM': 2.45,  # Curium
        'BK': 2.44,  # Berkelium
        'CF': 2.45,  # Californium
        'ES': 2.45,  # Einsteinium
        'FM': 2.45,  # Fermium
        'MD': 2.46,  # Mendelevium
        'NO': 2.46,  # Nobelium
        'LR': 2.46,  # Lawrencium
    }
    
    @classmethod
    def get_radius(cls, element, default_radius=1.5):
        """
        Get van der Waals radius for an element.
        
        Args:
            element: Element symbol (e.g., 'C', 'N', 'O')
            default_radius: Default radius if element not found (Angstroms)
            
        Returns:
            radius: van der Waals radius in Angstroms
        """
        # Convert to uppercase and handle common variations
        element_clean = str(element).upper().strip()
        
        return cls.VDW_RADII.get(element_clean, default_radius)
    
    @classmethod
    def assign_radii_to_atoms(cls, atoms, default_radius=1.5, verbose=True):
        """
        Assign van der Waals radii to MDAnalysis atoms based on their elements.
        Enhanced to handle common atom naming conventions.
        
        Args:
            atoms: MDAnalysis AtomGroup
            default_radius: Default radius for unknown elements (Angstroms)
            verbose: Print statistics about radius assignment
            
        Returns:
            radii: numpy array of radii in Angstroms
        """
        radii = np.zeros(len(atoms))
        element_counts = {}
        unknown_elements = set()
        
        # Define element aliases for common atom naming conventions
        element_aliases = {
            'CA+': 'CA', 'MG+': 'MG', 'NA+': 'NA', 'K+': 'K', 'CL-': 'CL',
            'SO4': 'S', 'PO4': 'P',
            # Common hydrogen variations
            'HA': 'H', 'HB': 'H', 'HG': 'H', 'HD': 'H', 'HE': 'H', 
            'HZ': 'H', 'HH': 'H', 'HN': 'H', 'H1': 'H', 'H2': 'H', 'H3': 'H',
            # Other common variations  
            'OG': 'O', 'OD': 'O', 'OE': 'O', 'OH': 'O',
            'NE': 'N', 'NH': 'N', 'NZ': 'N', 'ND': 'N',
            'SG': 'S', 'SD': 'S',
        }
        
        for i, atom in enumerate(atoms):
            # Try to get element from different possible attributes
            element = None
            
            # First try the element attribute if it exists
            if hasattr(atom, 'element') and atom.element and atom.element.strip():
                element = atom.element.strip()
            elif hasattr(atom, 'name') and atom.name:
                atom_name = atom.name.strip()
                # Try direct lookup in aliases first
                if atom_name in element_aliases:
                    element = element_aliases[atom_name]
                # Handle common patterns
                elif atom_name.startswith('H'):
                    element = 'H'
                elif atom_name.startswith('C'):
                    element = 'C'  
                elif atom_name.startswith('N'):
                    element = 'N'
                elif atom_name.startswith('O'):
                    element = 'O'
                elif atom_name.startswith('S'):
                    element = 'S'
                elif atom_name.startswith('P'):
                    element = 'P'
                else:
                    # Fallback: try first letter, then first two letters
                    element = atom_name[0].upper()
                    if element not in cls.VDW_RADII and len(atom_name) > 1:
                        element = atom_name[:2].upper()
                        
            elif hasattr(atom, 'type') and atom.type:
                atom_type = atom.type.strip()
                # Apply same logic to atom type
                if atom_type in element_aliases:
                    element = element_aliases[atom_type]
                elif atom_type.startswith('H'):
                    element = 'H'
                else:
                    element = atom_type[0].upper()
                    if element not in cls.VDW_RADII and len(atom_type) > 1:
                        element = atom_type[:2].upper()
            
            if element:
                radius = cls.get_radius(element, default_radius)
                if element.upper() not in cls.VDW_RADII:
                    unknown_elements.add(element.upper())
                    
                # Count elements for statistics
                element_key = element.upper()
                element_counts[element_key] = element_counts.get(element_key, 0) + 1
            else:
                radius = default_radius
                unknown_elements.add('UNKNOWN')
                element_counts['UNKNOWN'] = element_counts.get('UNKNOWN', 0) + 1
            
            radii[i] = radius
        
        if verbose:
            logger.info("Van der Waals radii assignment statistics:")
            for element, count in sorted(element_counts.items()):
                if element in cls.VDW_RADII:
                    radius = cls.VDW_RADII[element]
                    logger.info(f"  {element}: {count} atoms, radius = {radius:.2f} Å")
                else:
                    logger.info(f"  {element}: {count} atoms, radius = {default_radius:.2f} Å (default)")
            
            if unknown_elements:
                logger.warning(f"Unknown elements using default radius: {unknown_elements}")
        
        return radii


class PoreGeometry:
    """
    Class to handle different pore geometries.
    """
    
    @staticmethod
    def create_cylindrical_pore(X, Y, Z, pore_radius, membrane_half_thickness, 
                            corner_radius=None, chamfer_depth=None):
        """
        Create a cylindrical pore with optional edge chamfering.
        
        Args:
            X, Y, Z: Meshgrid coordinates (Å)
            pore_radius: Pore radius at center (Å)
            membrane_half_thickness: Half membrane thickness (Å)
            corner_radius: Additional radius at edges beyond pore_radius (Å). If None, creates regular cylindrical pore.
                        Edge radius = pore_radius + corner_radius
            chamfer_depth: Depth of chamfer cut from membrane edge (Å). If None, defaults to corner_radius
                        for 45-degree chamfer. If specified, creates custom chamfer depth.
        
        Returns:
            mask: Boolean array indicating membrane regions (True = membrane, False = pore)
        """
        import numpy as np
        
        R = np.sqrt(X**2 + Y**2)
        
        # Determine local pore radius based on Z position
        if corner_radius is None:
            # Regular cylindrical pore - use constant radius
            local_pore_radius = np.full_like(Z, pore_radius)
        else:
            # Calculate edge radius
            edge_radius = pore_radius + corner_radius
            
            if chamfer_depth is None:
                # Default: chamfer extends inward by corner_radius distance (45-degree slope)
                chamfer_depth = corner_radius
            
            # Linear chamfer only near edges
            local_pore_radius = np.full_like(Z, pore_radius)
            
            # Distance from membrane edges
            z_edge_dist = membrane_half_thickness - np.abs(Z)
            
            # Apply chamfer only within chamfer_depth from edges
            in_chamfer_zone = z_edge_dist < chamfer_depth
                
            # Linear interpolation in chamfer zone
            # At edge (z_edge_dist=0): radius = edge_radius
            # At chamfer_depth: radius = pore_radius
            chamfer_progress = z_edge_dist / chamfer_depth  # 0 at edge, 1 at chamfer_depth
            chamfer_progress = np.clip(chamfer_progress, 0, 1)
                
            # Apply chamfer
            local_pore_radius = np.where(in_chamfer_zone,
                                        edge_radius + (pore_radius - edge_radius) * chamfer_progress,
                                        pore_radius)
        
        # Create membrane mask
        membrane_mask = (np.abs(Z) < membrane_half_thickness) & (R > local_pore_radius)
        
        return membrane_mask
    
    @staticmethod
    def create_double_cone_pore(X, Y, Z, inner_radius, outer_radius, membrane_half_thickness):
        """
        Create a double cone (hourglass) pore.
        
        Args:
            X, Y, Z: Meshgrid coordinates (Å)
            inner_radius: Inner radius at membrane center (Å)
            outer_radius: Outer radius at membrane surfaces (Å)
            membrane_half_thickness: Half membrane thickness (Å)
            
        Returns:
            mask: Boolean array indicating membrane regions (True = membrane, False = pore)
        """
        R = np.sqrt(X**2 + Y**2)
        
        # Calculate radius as function of Z position
        # Linear interpolation from outer_radius at edges to inner_radius at center
        abs_z = np.abs(Z)
        
        # Only consider points within the membrane thickness
        in_membrane_z = abs_z < membrane_half_thickness
        
        # Calculate local pore radius based on z position
        # At z=0 (center): radius = inner_radius
        # At z=±membrane_half_thickness: radius = outer_radius
        z_fraction = abs_z / membrane_half_thickness
        local_pore_radius = inner_radius + (outer_radius - inner_radius) * z_fraction
        
        # Membrane exists where we're in the membrane region AND outside the local pore radius
        membrane_mask = in_membrane_z & (R > local_pore_radius)
        
        return membrane_mask


class VerticalMovementSEM:
    """
    SEM calculation for vertical movement of analyte through nanopore.
    Enhanced version with van der Waals radii support and multiple pore geometries.
    All units in Angstroms.
    
    This class:
    1. Creates a mathematical conductivity grid for membrane with pore (cylindrical or double cone)
    2. Loads a moving analyte structure
    3. Assigns van der Waals radii based on atomic elements
    4. Generates virtual trajectory by moving analyte vertically
    5. Calculates conductance at each position
    6. Solves FEM for current at each position
    """
    
    def __init__(self, 
                 moving_pdb,        # Centered analyte
                 pore_type="cylindrical",  # "cylindrical" or "double_cone"
                 pore_radius=100.0,  # Pore radius (Å) - for cylindrical or inner radius for double cone
                 outer_radius=None,  # Outer radius for double cone (Å) - if None, uses pore_radius * 1.5
                 corner_radius=0.0,  # Corner radius for cylindrical pore (Å)
                 membrane_thickness=200.0,  # Membrane thickness (Å)
                 z_start=150.0,     # Starting Z position (Å)
                 z_end=-150.0,      # Ending Z position (Å)
                 z_step=1.0,        # Step size (Å)
                 voltage=100.0,     # Applied voltage (mV)
                 output_prefix="vertical_movement",
                 box_dimensions=None,  # If None, auto-calculate
                 grid_resolution=1.0,  # Grid resolution (Å)
                 cutoff=4.1,          # Distance cutoff for analyte (Å)
                 bulk_conductivity=10.5,  # Bulk conductivity (S/m)
                 use_vdw_radii=True,  # Use van der Waals radii
                 default_radius=1.5): # Default radius for unknown elements (Å)
        
        self.moving_pdb = moving_pdb
        self.pore_type = pore_type.lower()
        self.pore_radius = pore_radius
        self.outer_radius = outer_radius if outer_radius is not None else pore_radius * 1.5
        self.corner_radius = corner_radius
        self.membrane_thickness = membrane_thickness
        self.z_start = z_start
        self.z_end = z_end
        self.z_step = z_step
        self.voltage = voltage / 1000.0  # Convert to V
        self.output_prefix = output_prefix
        self.grid_resolution = grid_resolution
        self.cutoff = cutoff
        self.bulk_conductivity = bulk_conductivity
        self.use_vdw_radii = use_vdw_radii
        self.default_radius = default_radius
        
        # Validate pore type
        if self.pore_type not in ["cylindrical", "double_cone"]:
            raise ValueError("pore_type must be 'cylindrical' or 'double_cone'")
        
        # For double cone, ensure outer_radius > pore_radius (inner_radius)
        if self.pore_type == "double_cone" and self.outer_radius <= self.pore_radius:
            raise ValueError("For double_cone pore, outer_radius must be greater than pore_radius (inner_radius)")
        
        # Load analyte structure
        logger.info("Loading analyte structure...")
        self.moving_universe = mda.Universe(moving_pdb)
        
        # Assign radii using van der Waals values or original method
        if self.use_vdw_radii:
            logger.info("Assigning van der Waals radii to atoms...")
            self.moving_radii = VanDerWaalsRadii.assign_radii_to_atoms(
                self.moving_universe.atoms, 
                default_radius=self.default_radius,
                verbose=True
            )
        else:
            # Use original method
            try:
                self.moving_radii = self.moving_universe.atoms.radii
                logger.info("Using radii from PDB file")
            except:
                self.moving_radii = np.ones(len(self.moving_universe.atoms)) * self.default_radius
                logger.info(f"Using default radius {self.default_radius} Å for all atoms")
        
        # Set or calculate box dimensions
        if box_dimensions is None:
            self.calculate_box_dimensions()
        else:
            self.box_dimensions = box_dimensions
            
        # Create base conductivity grid for membrane
        self.create_base_conductivity_grid()
        
        # Initialize conductivity model for analyte
        self.conductivity_model = SimpleConductivityModel(
            bulk_conductivity=bulk_conductivity,
            cutoff=cutoff  # Keep in Angstroms
        )
        
        # Setup FEniCS
        self.setup_fenics()
        
    def calculate_box_dimensions(self):
        """Auto-calculate box dimensions based on pore geometry."""
        # Use fixed dimensions that encompass the membrane and movement range
        padding = 20.0  # Å
        
        # XY dimensions based on pore size and padding
        if self.pore_type == "cylindrical":
            max_radius = self.pore_radius
        else:  # double_cone
            max_radius = self.outer_radius
            
        xy_size = max(150.0, max_radius * 3 + padding)
        
        # Z dimension based on movement range and membrane
        z_min = min(-self.membrane_thickness/2 - padding, self.z_end - padding)
        z_max = max(self.membrane_thickness/2 + padding, self.z_start + padding)
        
        self.box_dimensions = {
            'x': (-xy_size, xy_size),
            'y': (-xy_size, xy_size),
            'z': (z_min, z_max)
        }
        
        logger.info(f"Box dimensions: X={self.box_dimensions['x']}, "
                   f"Y={self.box_dimensions['y']}, Z={self.box_dimensions['z']}")
    
    def create_base_conductivity_grid(self):
        """
        Create a mathematical conductivity grid for membrane with specified pore geometry.
        This replaces loading a membrane PDB file.
        """
        logger.info(f"Creating base conductivity grid for {self.pore_type} pore...")
        
        # Create grid in Angstroms
        x_range = np.arange(self.box_dimensions['x'][0], 
                           self.box_dimensions['x'][1], 
                           self.grid_resolution)
        y_range = np.arange(self.box_dimensions['y'][0], 
                           self.box_dimensions['y'][1], 
                           self.grid_resolution)
        z_range = np.arange(self.box_dimensions['z'][0], 
                           self.box_dimensions['z'][1], 
                           self.grid_resolution)
        
        # Create meshgrid
        X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing='ij')
        
        # Initialize with bulk conductivity
        self.base_conductivity = np.ones(X.shape) * self.bulk_conductivity
        
        # Define membrane region based on pore type
        membrane_half_thickness = self.membrane_thickness / 2  # Keep in Angstroms
        
        if self.pore_type == "cylindrical":
            pore_radius = self.pore_radius  # Keep in Angstroms
            corner_radius = self.corner_radius  # Keep in Angstroms
            
            membrane_mask = PoreGeometry.create_cylindrical_pore(
                X, Y, Z, pore_radius, membrane_half_thickness, corner_radius
            )
            
            logger.info(f"Cylindrical pore - radius: {self.pore_radius} Å, "
                       f"corner radius: {self.corner_radius} Å")
            
        elif self.pore_type == "double_cone":
            inner_radius = self.pore_radius  # Keep in Angstroms
            outer_radius = self.outer_radius  # Keep in Angstroms
            
            membrane_mask = PoreGeometry.create_double_cone_pore(
                X, Y, Z, inner_radius, outer_radius, membrane_half_thickness
            )
            
            logger.info(f"Double cone pore - inner radius: {self.pore_radius} Å, "
                       f"outer radius: {self.outer_radius} Å")
        
        # Set low conductivity in membrane regions
        self.base_conductivity[membrane_mask] = 0.0001  # Very low conductivity
        
        # Store grid edges for interpolation
        self.base_conductivity_edges = (x_range, y_range, z_range)
        
        # Create interpolator for base conductivity
        self.base_cond_interp = RegularGridInterpolator(
            self.base_conductivity_edges,
            self.base_conductivity,
            bounds_error=False,
            fill_value=self.bulk_conductivity
        )
        
        logger.info(f"Created base conductivity grid with shape {self.base_conductivity.shape}")
        logger.info(f"Membrane thickness: {self.membrane_thickness} Å ({membrane_half_thickness*2} Å)")
        
        # Log some statistics
        n_membrane = np.sum(membrane_mask)
        n_total = X.size
        logger.info(f"Membrane occupies {n_membrane}/{n_total} grid points "
                   f"({100*n_membrane/n_total:.1f}%)")
    
    def setup_fenics(self):
        """Initialize FEniCS solver components."""
        # Create mesh
        logger.info("Creating FEniCS mesh...")
        
        nx = int((self.box_dimensions['x'][1] - self.box_dimensions['x'][0]) / self.grid_resolution)
        ny = int((self.box_dimensions['y'][1] - self.box_dimensions['y'][0]) / self.grid_resolution)
        nz = int((self.box_dimensions['z'][1] - self.box_dimensions['z'][0]) / self.grid_resolution)
        
        # Ensure reasonable mesh size
        max_elements = 400
        if max(nx, ny, nz) > max_elements:
            scale_factor = max_elements / max(nx, ny, nz)
            nx = int(nx * scale_factor)
            ny = int(ny * scale_factor)
            nz = int(nz * scale_factor)
            logger.warning(f"Reduced mesh resolution to {nx}x{ny}x{nz} for computational efficiency")
        
        self.mesh = BoxMesh(
            Point(self.box_dimensions['x'][0], self.box_dimensions['y'][0], self.box_dimensions['z'][0]),
            Point(self.box_dimensions['x'][1], self.box_dimensions['y'][1], self.box_dimensions['z'][1]),
            nx, ny, nz
        )
        
        logger.info(f"Created mesh with {nx}x{ny}x{nz} elements")
        
        # Define function spaces
        self.V = FunctionSpace(self.mesh, 'P', 1)  # For potential
        self.F = FunctionSpace(self.mesh, 'CG', 1)  # For conductivity
        
        # Define boundary conditions
        z_min = self.box_dimensions['z'][0]
        z_max = self.box_dimensions['z'][1]
        
        ground = CompiledSubDomain('on_boundary && near(x[2], z_min, tol)', 
                                  z_min=z_min, tol=1e-6)
        terminal = CompiledSubDomain('on_boundary && near(x[2], z_max, tol)', 
                                    z_max=z_max, tol=1e-6)
        
        self.bc_ground = DirichletBC(self.V, Constant(0), ground)
        self.bc_terminal = DirichletBC(self.V, Constant(self.voltage), terminal)
        self.bcs = [self.bc_ground, self.bc_terminal]
        
        # Setup boundary markers for flux calculation
        boundary_parts = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        ground.mark(boundary_parts, 1)
        terminal.mark(boundary_parts, 2)
        self.ds = Measure('ds', domain=self.mesh, subdomain_data=boundary_parts)
        
        # Define variational problem
        u = TrialFunction(self.V)
        v = TestFunction(self.V)
        self.sig = Function(self.F)
        self.u1 = Function(self.V)
        
        f = Constant(0)
        DE = self.sig * dot(grad(u), grad(v)) * dx - f * v * dx
        self.a, self.L = lhs(DE), rhs(DE)
        
        # Setup solver
        self.solver = KrylovSolver("gmres", "amg")
        self.solver.parameters["relative_tolerance"] = 1e-8
        self.solver.parameters["maximum_iterations"] = 20000
        
        logger.info("FEniCS setup complete")
    
    def get_conductivity_at_position(self, z_position):
        """
        Calculate conductivity field when moving atoms are at given z position.
        Enhanced to use van der Waals radii.
        
        Args:
            z_position: Z coordinate to place the center of mass of moving atoms
            
        Returns:
            Conductivity field as numpy array
        """
        # Get mesh coordinates
        dofmap = self.F.dofmap()
        my_first, my_last = dofmap.ownership_range()
        n = self.F.dim()
        d = self.mesh.geometry().dim()
        F_dof_coordinates = self.F.tabulate_dof_coordinates()
        F_dof_coordinates.resize((n, d))
        
        # Mesh coordinates are already in Angstroms
        mesh_coords = F_dof_coordinates
        
        # Get base conductivity at mesh points
        base_cond = self.base_cond_interp(mesh_coords)
        
        # Handle NaN values (outside interpolation domain)
        z_coords = mesh_coords[:, 2]
        membrane_half_thick = self.membrane_thickness / 2  # Angstroms
        
        # Points outside interpolation should have bulk conductivity
        nan_mask = np.isnan(base_cond)
        base_cond[nan_mask] = self.bulk_conductivity
        
        # Now modify conductivity based on analyte position
        moving_atoms = self.moving_universe.atoms
        
        # Calculate current COM of moving atoms
        moving_com = moving_atoms.center_of_mass()
        displacement = np.array([0, 0, z_position - moving_com[2]])
        
        # Get displaced positions in Angstroms
        moving_positions = moving_atoms.positions + displacement
        
        # Get radii in Angstroms (use the pre-assigned VdW radii)
        moving_radii = self.moving_radii
        
        # Calculate modification due to analyte
        analyte_cond = self.calculate_analyte_conductivity_modification(
            mesh_coords, moving_positions, moving_radii, base_cond
        )
        
        return analyte_cond
    
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
            distance_upper_bound=cutoff + np.max(atom_radii)
        )
        
        # Adjust distances for atom radii
        valid_mask = ~np.isinf(distances)
        if np.any(valid_mask):
            distances[valid_mask] -= atom_radii[indices[valid_mask]]
            distances[valid_mask] = np.maximum(distances[valid_mask], 0)
            
            # Apply conductivity model to get modulation factor
            analyte_cond = self.conductivity_model(distances[valid_mask])
            
            # Modify conductivity: use minimum of base and analyte-modified conductivity
            conductivity[valid_mask] = np.minimum(
                conductivity[valid_mask],
                analyte_cond
            )
        
        return conductivity
    
    def solve_for_current(self, conductivity):
        """
        Solve FEM problem for given conductivity field and return current.
        
        Args:
            conductivity: Conductivity field values at mesh DOFs
            
        Returns:
            current: Calculated current (A)
        """
        # Load conductivity into FEniCS function
        vec = self.sig.vector()
        vec.set_local(conductivity)
        vec.apply('insert')
        
        # Assemble and solve
        A, b = assemble_system(self.a, self.L, self.bcs)
        self.solver.set_operator(A)
        self.solver.solve(self.u1.vector(), b)
        
        # Calculate flux at boundaries
        flux_in = dot(Constant((0, 0, 1)), self.sig * nabla_grad(self.u1)) * self.ds(1)
        flux_out = dot(Constant((0, 0, 1)), self.sig * nabla_grad(self.u1)) * self.ds(2)
        
        current_in = assemble(flux_in)
        current_out = assemble(flux_out)
        
        # Return average magnitude
        return abs(current_out)
    
    def calculate_open_pore_current(self):
        """Calculate baseline current with no analyte in the system."""
        logger.info("Calculating open pore current...")
        
        # Get mesh coordinates
        dofmap = self.F.dofmap()
        n = self.F.dim()
        d = self.mesh.geometry().dim()
        F_dof_coordinates = self.F.tabulate_dof_coordinates()
        F_dof_coordinates.resize((n, d))
        mesh_coords = F_dof_coordinates  # Already in Angstroms
        
        # Get base conductivity
        base_cond = self.base_cond_interp(mesh_coords)
        nan_mask = np.isnan(base_cond)
        base_cond[nan_mask] = self.bulk_conductivity
        
        # Solve for current
        open_current = self.solve_for_current(base_cond)
        logger.info(f"Open pore current: {open_current:.6e} A")
        
        return open_current
    
    def run(self):
        """
        Run the complete vertical movement simulation.
        
        Returns:
            results: Dictionary with z_positions, currents, and normalized currents
        """
        # First calculate open pore current
        open_current = self.calculate_open_pore_current()
        
        # Calculate number of steps
        num_steps = int(abs(self.z_end - self.z_start) / self.z_step) + 1
        z_positions = np.linspace(self.z_start, self.z_end, num_steps)
        
        logger.info(f"Running simulation with {num_steps} positions")
        logger.info(f"Z range: {self.z_start} to {self.z_end} Å, step: {self.z_step} Å")
        logger.info(f"Pore type: {self.pore_type}")
        if self.pore_type == "cylindrical":
            logger.info(f"Pore radius: {self.pore_radius} Å")
            if self.corner_radius > 0:
                logger.info(f"Corner radius: {self.corner_radius} Å")
        else:  # double_cone
            logger.info(f"Inner radius: {self.pore_radius} Å")
            logger.info(f"Outer radius: {self.outer_radius} Å")
        
        if self.use_vdw_radii:
            logger.info("Using van der Waals radii for accurate atomic volumes")
        else:
            logger.info(f"Using uniform radius of {self.default_radius} Å for all atoms")
        
        currents = []
        
        for i, z_pos in enumerate(z_positions):
            logger.info(f"Processing position {i+1}/{num_steps}: Z = {z_pos:.1f} Å")
            
            # Calculate conductivity field
            conductivity = self.get_conductivity_at_position(z_pos)
            
            # Solve for current
            current = self.solve_for_current(conductivity)
            currents.append(current)
            
            # Calculate blockage
            blockage = (1 - current/open_current) * 100
            logger.info(f"  Current: {current:.6e} A (blockage: {blockage:.1f}%)")
            
            # Save intermediate result
            with open(f"{self.output_prefix}_position_{i:04d}.dat", 'w') as f:
                f.write(f"{z_pos} {current} {blockage}\n")
        
        currents = np.array(currents)
        normalized_currents = currents / open_current
        blockages = (1 - normalized_currents) * 100
        
        # Save final results
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
            'use_vdw_radii': self.use_vdw_radii
        }
        
        header = f"Z_position(Angstrom) Current(Ampere) Normalized_Current Blockage(%) # Pore_type: {self.pore_type}"
        if self.use_vdw_radii:
            header += " - Using van der Waals radii"
        
        np.savetxt(f"{self.output_prefix}_results.txt", 
                   np.column_stack([z_positions, currents, normalized_currents, blockages]),
                   header=header)
        
        logger.info("Simulation complete!")
        logger.info(f"Maximum blockage: {np.max(blockages):.1f}% at Z = {z_positions[np.argmax(blockages)]:.1f} Å")
        
        return results

    def get_conductivity_grid_for_preview(self, z_position):
        """
        Get 2D conductivity slice for visualization at given z position.
        Enhanced to work with different pore geometries.
        
        Args:
            z_position: Z coordinate to place the center of mass of moving atoms
            
        Returns:
            x_coords, z_coords, conductivity_2d: Arrays for plotting
        """
        # Create a 2D slice through the middle of the Y dimension for visualization
        y_middle = (self.box_dimensions['y'][0] + self.box_dimensions['y'][1]) / 2
        
        # Create coordinate arrays for the slice
        x_coords = np.arange(self.box_dimensions['x'][0], 
                            self.box_dimensions['x'][1], 
                            self.grid_resolution)
        z_coords = np.arange(self.box_dimensions['z'][0], 
                            self.box_dimensions['z'][1], 
                            self.grid_resolution)
        
        # Create meshgrid for the slice
        X_slice, Z_slice = np.meshgrid(x_coords, z_coords, indexing='ij')
        Y_slice = np.full_like(X_slice, y_middle)
        
        # Create coordinate array for interpolation (in Angstroms)
        slice_coords = np.column_stack([
            X_slice.ravel(),  # Keep in Angstroms
            Y_slice.ravel(),  # Keep in Angstroms
            Z_slice.ravel()   # Keep in Angstroms
        ])
        
        # Get base conductivity at slice points
        base_cond = self.base_cond_interp(slice_coords)
        
        # Handle NaN values
        nan_mask = np.isnan(base_cond)
        base_cond[nan_mask] = self.bulk_conductivity
        
        # Now modify conductivity based on analyte position
        moving_atoms = self.moving_universe.atoms
        
        # Calculate current COM of moving atoms
        moving_com = moving_atoms.center_of_mass()
        displacement = np.array([0, 0, z_position - moving_com[2]])
        
        # Get displaced positions in Angstroms
        moving_positions = moving_atoms.positions + displacement
        
        # Get radii in Angstroms (use the pre-assigned VdW radii)
        moving_radii = self.moving_radii
        
        # Calculate modification due to analyte
        analyte_cond = self.calculate_analyte_conductivity_modification(
            slice_coords, moving_positions, moving_radii, base_cond
        )
        
        # Reshape back to 2D
        conductivity_2d = analyte_cond.reshape(X_slice.shape)
        
        return x_coords, z_coords, conductivity_2d

    def preview_analyte_movement(self, num_preview_frames=4, save_plots=True):
        """
        Create preview plots showing the analyte movement through the nanopore.
        Enhanced to show different pore geometries.
        
        Args:
            num_preview_frames: Number of frames to show
            save_plots: Whether to save the plots
        """
        logger.info(f"Creating {num_preview_frames} preview frames of analyte movement...")
        logger.info(f"Pore type: {self.pore_type}")
        
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("Matplotlib not available, cannot create preview plots")
            return
        
        # Calculate total number of steps based on z_start, z_end, and z_step
        total_steps = int(abs(self.z_end - self.z_start) / self.z_step) + 1
        z_positions = np.linspace(self.z_start, self.z_end, total_steps)
        
        logger.info(f"Total available steps: {total_steps} (from Z={self.z_start} to Z={self.z_end} with step={self.z_step})")
        logger.info(f"Selecting {num_preview_frames} frames from {total_steps} total steps")
        
        # Select frames to preview - evenly distribute across all positions
        if num_preview_frames >= total_steps:
            # If requesting more frames than available positions, use all positions
            frame_indices = list(range(total_steps))
            logger.info(f"Using all {total_steps} available positions (requested {num_preview_frames})")
        else:
            # Evenly distribute the requested number of frames across the total range
            frame_indices = np.linspace(0, total_steps-1, num_preview_frames, dtype=int)
            frame_indices = list(frame_indices)  # Convert to list for consistency
            
            # Log which specific positions are being used
            selected_z_positions = [z_positions[idx] for idx in frame_indices]
            logger.info(f"Selected Z positions: {[f'{z:.1f}' for z in selected_z_positions]}")
        
        for i, frame in enumerate(frame_indices):
            z_pos = z_positions[frame]
            logger.info(f"Creating preview frame {i+1}/{len(frame_indices)} at Z = {z_pos:.1f} Å")
            
            # Get conductivity grid for this position
            x_coords, z_coords, cond = self.get_conductivity_grid_for_preview(z_pos)
            
            # Create the plot
            plt.figure(figsize=(12, 8))
            im = plt.imshow(cond.T, extent=[z_coords[0], z_coords[-1], x_coords[0], x_coords[-1]], 
                        aspect='auto', origin='lower', cmap='viridis')
            plt.colorbar(im, label='Conductivity (S/m)')
            plt.xlabel('X position (Å)')
            plt.ylabel('Z position (Å)')
            
            # Create title with pore geometry information
            if self.pore_type == "cylindrical":
                if self.corner_radius > 0:
                    title = f'Conductivity Side View - {self.pore_type.title()} Pore (R={self.pore_radius:.0f}Å, Corner R={self.corner_radius:.0f}Å)\nAnalyte at Z = {z_pos:.1f} Å'
                else:
                    title = f'Conductivity Side View - {self.pore_type.title()} Pore (R={self.pore_radius:.0f}Å)\nAnalyte at Z = {z_pos:.1f} Å'
            else:  # double_cone
                title = f'Conductivity Side View - {self.pore_type.title()} Pore (Inner R={self.pore_radius:.0f}Å, Outer R={self.outer_radius:.0f}Å)\nAnalyte at Z = {z_pos:.1f} Å'
            
            plt.title(title)
            
            # Add membrane boundaries
            membrane_half_thick = self.membrane_thickness / 2
            plt.axhline(y=membrane_half_thick, color='red', linestyle='--', alpha=0.7, label='Membrane')
            plt.axhline(y=-membrane_half_thick, color='red', linestyle='--', alpha=0.7)
            
            # Add pore boundaries based on geometry type
            if self.pore_type == "cylindrical":
                # Simple cylindrical pore boundaries
                plt.axvline(x=self.pore_radius, color='red', linestyle=':', alpha=0.5, 
                        label=f'Pore (R={self.pore_radius:.0f}Å)')
                plt.axvline(x=-self.pore_radius, color='red', linestyle=':', alpha=0.5)
                    
            else:  # double_cone
                # Extract the actual boundary from the conductivity data
                z_membrane_coords = z_coords[(z_coords >= -membrane_half_thick) & (z_coords <= membrane_half_thick)]
                if len(z_membrane_coords) > 0:
                    # Also plot the theoretical boundary for comparison
                    abs_z = np.abs(z_membrane_coords)
                    z_fraction = abs_z / membrane_half_thick
                    local_pore_radius = self.pore_radius + (self.outer_radius - self.pore_radius) * z_fraction
                    
                    plt.plot(local_pore_radius, z_membrane_coords, 'red', linestyle=':', alpha=0.7, 
                            linewidth=2, label=f'Theoretical boundary')
                    plt.plot(-local_pore_radius, z_membrane_coords, 'red', linestyle=':', alpha=0.7, 
                            linewidth=2)
                    
                    # Show inner and outer radii with horizontal lines at membrane edges
                    plt.axvline(x=self.outer_radius, color='blue', linestyle='-.', alpha=0.5, 
                            label=f'Outer R={self.outer_radius:.0f}Å')
                    plt.axvline(x=-self.outer_radius, color='blue', linestyle='-.', alpha=0.5)
                    plt.axvline(x=self.pore_radius, color='green', linestyle='-.', alpha=0.5, 
                            label=f'Inner R={self.pore_radius:.0f}Å')
                    plt.axvline(x=-self.pore_radius, color='green', linestyle='-.', alpha=0.5)
                        
            # Add analyte position indicator
            plt.axhline(y=z_pos, color='white', linestyle='-', alpha=0.8, linewidth=2, 
                    label=f'Analyte COM')
            
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.grid(True, alpha=0.3)
            
            if save_plots:
                filename = f'{self.output_prefix}_{self.pore_type}_preview_frame_{frame:04d}_z_{z_pos:.1f}A.png'
                plt.savefig(filename, dpi=150, bbox_inches='tight')
                logger.info(f"Saved preview frame: {filename}")
                
                # Also save data
                data_filename = f'{self.output_prefix}_{self.pore_type}_preview_frame_{frame:04d}_z_{z_pos:.1f}A.dat'
                np.savetxt(data_filename, cond)
            
            plt.show()
            plt.close()
        
        logger.info("Preview frames created successfully!")
        
    def preview_only(self, num_preview_frames=4):
        """
        Create only preview plots without running SEM calculations.
        
        Args:
            num_preview_frames: Number of preview frames to create
        """
        logger.info("Running preview-only mode (no SEM calculations)")
        logger.info(f"Pore geometry: {self.pore_type}")
        if self.pore_type == "cylindrical":
            logger.info(f"Pore radius: {self.pore_radius} Å")
            if self.corner_radius > 0:
                logger.info(f"Corner radius: {self.corner_radius} Å")
        else:  # double_cone
            logger.info(f"Inner radius: {self.pore_radius} Å")
            logger.info(f"Outer radius: {self.outer_radius} Å")
        
        self.preview_analyte_movement(num_preview_frames=num_preview_frames)


class SimpleConductivityModel:
    """
    Simple conductivity model based on distance from atoms.
    All units in Angstroms.
    """
    
    def __init__(self, bulk_conductivity=1.12, cutoff=4.1):
        """
        Args:
            bulk_conductivity: Bulk conductivity (S/m)
            cutoff: Distance cutoff (Å)
        """
        self.bulk_conductivity = bulk_conductivity
        self.cutoff = cutoff
        self.min_distance = 1.3  # Å
        
    def __call__(self, distances):
        """
        Calculate conductivity based on distance from atoms.
        
        Args:
            distances: Array of distances (Å)
            
        Returns:
            conductivity: Array of conductivity values (S/m)
        """
        # Linear interpolation between min_distance and cutoff
        # 0 conductivity at min_distance, bulk conductivity at cutoff
        
        conductivity = np.zeros_like(distances)
        min_conductivity = 0.0000001 * self.bulk_conductivity
        # Inside cutoff
        mask = distances < self.cutoff
        if np.any(mask):
            # Linear interpolation
            fraction = (distances[mask] - self.min_distance) / (self.cutoff - self.min_distance)
            fraction = np.clip(fraction, 0, 1)
            conductivity[mask] = min_conductivity + fraction * self.bulk_conductivity
        
        # Outside cutoff
        conductivity[distances >= self.cutoff] = self.bulk_conductivity
        
        return conductivity


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
        logger.info(f"Loaded configuration from {config_path}")
        return config
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing JSON configuration: {e}")
        sys.exit(1)


def create_sem_from_config(config):
    """
    Create VerticalMovementSEM instance from configuration dictionary.
    
    Args:
        config: Configuration dictionary loaded from JSON
        
    Returns:
        sem: VerticalMovementSEM instance
        config: Configuration dictionary (returned for convenience)
    """
    logger.info("Creating SEM instance from configuration...")
    
    # Extract parameters from config
    moving_pdb = config["input"]["moving_pdb"]
    
    # Pore geometry parameters
    pore_geom = config["pore_geometry"]
    pore_type = pore_geom["pore_type"].lower()  # Convert to lowercase
    pore_radius = pore_geom["pore_radius"]
    corner_radius = pore_geom.get("corner_radius", 0.0)
    membrane_thickness = pore_geom["membrane_thickness"]
    
    # For double cone, get outer_radius if specified
    outer_radius = pore_geom.get("outer_radius", None)
    
    # Simulation parameters
    sim = config["simulation"]
    voltage = sim["voltage"]
    bulk_conductivity = sim["bulk_conductivity"]
    grid_resolution = sim["grid_resolution"]
    use_vdw_radii = sim["use_vdw_radii"]
    default_radius = sim["default_radius"]
    
    # Movement parameters
    movement = config["movement"]
    z_start = movement["z_start"]
    z_end = movement["z_end"]
    z_step = movement["z_step"]
    
    # Output parameters
    output = config["output"]
    output_prefix = output["output_prefix"]
    
    # Box dimensions (optional)
    box_dims = config.get("box_dimensions", None)
    if box_dims:
        box_dimensions = {
            'x': tuple(box_dims["x"]),
            'y': tuple(box_dims["y"]),
            'z': tuple(box_dims["z"])
        }
        logger.info("Using box dimensions from configuration")
    else:
        box_dimensions = None
        logger.info("Box dimensions will be auto-calculated")
    
    # Create SEM instance
    sem = VerticalMovementSEM(
        moving_pdb=moving_pdb,
        pore_type=pore_type,
        pore_radius=pore_radius,
        outer_radius=outer_radius,
        corner_radius=corner_radius,
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
        default_radius=default_radius
    )
    
    logger.info("SEM instance created successfully")
    return sem, config


def print_config_summary(config):
    """
    Print a summary of the loaded configuration.
    
    Args:
        config: Configuration dictionary
    """
    logger.info("Configuration Summary:")
    logger.info(f"  Input PDB: {config['input']['moving_pdb']}")
    
    pore_geom = config["pore_geometry"]
    logger.info(f"  Pore Type: {pore_geom['pore_type']}")
    logger.info(f"  Pore Radius: {pore_geom['pore_radius']} Å")
    if pore_geom.get("corner_radius", 0) > 0:
        logger.info(f"  Corner Radius: {pore_geom['corner_radius']} Å")
    if "outer_radius" in pore_geom:
        logger.info(f"  Outer Radius: {pore_geom['outer_radius']} Å")
    logger.info(f"  Membrane Thickness: {pore_geom['membrane_thickness']} Å")
    
    sim = config["simulation"]
    logger.info(f"  Voltage: {sim['voltage']} mV")
    logger.info(f"  Bulk Conductivity: {sim['bulk_conductivity']} S/m")
    logger.info(f"  Use VdW Radii: {sim['use_vdw_radii']}")
    
    movement = config["movement"]
    logger.info(f"  Z Range: {movement['z_start']} to {movement['z_end']} Å")
    logger.info(f"  Z Step: {movement['z_step']} Å")
    
    output = config["output"]
    logger.info(f"  Output Prefix: {output['output_prefix']}")
    logger.info(f"  Preview Frames: {output.get('preview_frames', 0)}")


def main():
    """
    Main function to run SEM with JSON configuration.
    Supports two modes: 'run' and 'preview_only'
    """
    parser = argparse.ArgumentParser(
        description='Run SEM calculation with JSON configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sem_script.py config.json run           # Run full simulation
  python sem_script.py config.json preview_only  # Generate preview plots only
        """
    )
    parser.add_argument('config', help='Path to JSON configuration file')
    parser.add_argument('mode', choices=['run', 'preview_only'], 
                       help='Run mode: "run" for full simulation, "preview_only" for visualization')
    
    args = parser.parse_args()
    
    try:
        # Load configuration
        config = load_config(args.config)
        
        # Print configuration summary
        print_config_summary(config)
        
        # Create SEM instance
        sem, config = create_sem_from_config(config)
        
        # Run based on mode
        if args.mode == 'preview_only':
            logger.info("Running in preview-only mode")
            preview_frames = config["output"].get("preview_frames", 4)
            sem.preview_only(num_preview_frames=preview_frames)
            
        elif args.mode == 'run':
            logger.info("Running full simulation")
            # Run full simulation
            results = sem.run()
            
            # Optionally create preview frames after simulation
            # preview_frames = config["output"].get("preview_frames", 0)
            # if preview_frames > 0:
            #     logger.info(f"Creating {preview_frames} preview frames...")
            #     sem.preview_analyte_movement(num_preview_frames=preview_frames)
        
        logger.info("Execution completed successfully!")
        
    except Exception as e:
        logger.error(f"Error during execution: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()