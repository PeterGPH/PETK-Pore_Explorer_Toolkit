#!/usr/bin/env python3
"""
PDB to XYZ converter script for SEM package.
Python equivalent of VMD TCL script for PDB coordinate extraction with van der Waals radii support.

Usage: python pdb2xyz.py pore.pdb [output.xyz]
Example: python pdb2xyz.py protein.pdb protein.xyz

Can be used as standalone script or as part of the SEM package.
"""

# Package compatibility - allows script to work both as module and standalone
if __name__ == "__main__" and __package__ is None:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent))

import sys
import numpy as np
import logging

# Set up logging to suppress some MDAnalysis warnings but keep important info
logging.basicConfig(level=logging.WARNING)

try:
    import MDAnalysis as mda
    from MDAnalysis.core.selection import SelectionError
except ImportError:
    print("Error: MDAnalysis not found. Install with: pip install MDAnalysis")
    sys.exit(1)

# Import VanDerWaalsRadii from the package
try:
    from van_der_waals import VanDerWaalsRadii
except ImportError:
    try:
        from sem.van_der_waals import VanDerWaalsRadii
    except ImportError:
        try:
            from ..van_der_waals import VanDerWaalsRadii
        except ImportError:
            print("Error: Could not import VanDerWaalsRadii class")
            print("Make sure van_der_waals.py is available or install the sem package")
            sys.exit(1)

def extract_pore_coordinates(pdb_file, output_file="pore.xyz"):
    """
    Extract coordinates and radius from PDB file, excluding water and ions
    
    Args:
        pdb_file (str): Path to input PDB file
        output_file (str): Path to output XYZ file
    """
    
    try:
        # Initialize PDB parser
        print(f"Loading PDB file: {pdb_file}")
        
        # Parse PDB file
        u = mda.Universe(pdb_file)
        
        # Select atoms excluding water and ions
        # Common water residue names: HOH, WAT, TIP3, SOL
        # Common ion names: NA, CL, K, MG, CA, ZN, etc.
        selection_string = (
            "not (resname HOH or resname WAT or resname TIP3 or resname SOL or "
            "resname NA or resname CL or resname K or resname MG or resname CA or "
            "resname ZN or resname FE or resname BR or resname I)"
        )
        
        try:
            sel_pore = u.select_atoms(selection_string)
        except SelectionError:
            print("Warning: Selection failed, trying simpler selection...")
            # Fallback to protein selection if the above fails
            sel_pore = u.select_atoms("protein")
        
        num_atoms = len(sel_pore)
        print(f"Pore atoms: {num_atoms}")
        
        if num_atoms == 0:
            print("Warning: No atoms selected. Check your PDB file.")
            return
        
        # Get coordinates
        coordinates = sel_pore.positions
        
        # Get radius information using VanDerWaalsRadii
        radius = get_radius_info(sel_pore)
        
        # Write coordinates to output file
        write_output(coordinates, radius, output_file)
        print(f"Pore coordinates written to {output_file}")
        
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        sys.exit(1)

def get_radius_info(selection):
    """
    Attempt to get radius information for atoms
    
    Args:
        selection: MDAnalysis atom selection
        
    Returns:
        numpy array of radius values
    """
    
    # Try to get radius from occupancy field first (like VMD reference)
    try:
        if hasattr(selection, 'occupancies'):
            occupancy = selection.occupancies
            # Check if occupancies look like radius values (positive, reasonable range)
            # and have some variation (not all the same value)
            if (np.all(occupancy > 0) and 
                np.all(occupancy < 10) and 
                len(np.unique(occupancy)) > 1 and
                np.mean(occupancy) > 0.5):  # Most radii should be > 0.5 Å
                print("Using occupancy column as radius")
                return occupancy
    except:
        pass
    
    # Try tempfactor as backup
    try:
        if hasattr(selection, 'tempfactors'):
            radius = selection.tempfactors
            # Check if tempfactors look like radius values (positive, reasonable range)
            if (np.all(radius > 0) and 
                np.all(radius < 10) and 
                len(np.unique(radius)) > 1 and
                np.mean(radius) > 0.5):
                print("Using tempfactor column as radius")
                return radius
    except:
        pass
    
    # Fallback: use van der Waals radii based on element
    print("Using van der Waals radii based on element types")
    radii = VanDerWaalsRadii.assign_radii_to_atoms(selection, default_radius=1.70, verbose=False)
    return radii

def write_output(coordinates, radius, output_file):
    """
    Write coordinates and radius to output file
    
    Args:
        coordinates: numpy array of xyz coordinates
        radius: numpy array of radius values
        output_file: output file path
    """
    
    with open(output_file, 'w') as f:
        for i, (xyz, r) in enumerate(zip(coordinates, radius)):
            f.write(f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f} {r:.3f}\n")

def main():
    """Main function to handle command line arguments"""
    
    if len(sys.argv) < 2:
        print("Usage: python pdb2xyz.py pore.pdb [output.xyz]")
        print("Example: python pdb2xyz.py protein.pdb protein.xyz")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    
    # Optional: allow custom output file as second argument
    output_file = sys.argv[2] if len(sys.argv) > 2 else "pore.xyz"
    
    print(f"Pore PDB: {pdb_file}")
    extract_pore_coordinates(pdb_file, output_file)

if __name__ == "__main__":
    main()