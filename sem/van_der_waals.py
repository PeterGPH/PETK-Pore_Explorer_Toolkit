"""
Van der Waals radii database and assignment functions.
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)

class VanDerWaalsRadii:
    """
    Van der Waals radii database for common elements in Angstroms.
    Based on Bondi (1964) and more recent literature values.
    """

    # Atom-name specific overrides for common biomolecular conventions
    ATOM_NAME_RADII = {
        "HN": 1.00,
        "HA1": 1.00,
        "HA2": 1.00,
        "CA": 1.60,  # Alpha carbon
    }

    # Standard van der Waals radii in Angstroms
    VDW_RADII = {
        'H': 1.00,   # Hydrogen
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
    def guess_element_from_pdb_line(cls, line: str) -> str:
        """
        Guess an element symbol from a PDB record line.
        """
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element:
            return element

        name_field = line[12:16].strip() if len(line) >= 16 else ""
        if not name_field:
            return ""

        token = "".join(ch for ch in name_field if not ch.isdigit()).upper()
        if not token:
            token = name_field.strip().upper()

        if len(token) >= 2 and token[:2] in cls.VDW_RADII:
            return token[:2]
        if token[:1] in cls.VDW_RADII:
            return token[:1]
        return token[:1]
    
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
                atom_name_key = atom.name.strip().upper() if hasattr(atom, 'name') and atom.name else ""
                element_key = element.upper()
                if (
                    atom_name_key
                    and atom_name_key in cls.ATOM_NAME_RADII
                    and not (element_key == atom_name_key and element_key in cls.VDW_RADII)
                ):
                    radius = cls.ATOM_NAME_RADII[atom_name_key]
                else:
                    radius = cls.get_radius(element_key, default_radius)
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
