"""
Conductivity models for SEM calculations.
Enhanced to include charge effects using empirical Debye-Hückel theory.
"""

import numpy as np
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

class SimpleConductivityModel:
    """
    Simple conductivity model based on distance from atoms.
    All units in Angstroms.
    """
    
    def __init__(self, bulk_conductivity=11.2, cutoff=4.1):
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
            conductivity: Array of conductivity values (S/ms)
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

class ChargeAwareConductivityModel(SimpleConductivityModel):
    """
    Charge-aware extension: Modulates steric conductivity with cosh(beta * phi) for ion density effects.
    """
    def __init__(self, bulk_conductivity=10.5, cutoff=4.1, charge_clip=2.0):
        super().__init__(bulk_conductivity=bulk_conductivity, cutoff=cutoff)
        self.charge_clip = charge_clip
    
    def __call__(self, distances, potentials):
        # Get base steric conductivity
        steric_cond = super().__call__(distances)
        
        # Charge modulation (cosh for symmetric electrolyte enhancement)
        charge_factor = np.cosh(potentials)  # potentials are dimensionless beta * phi
        charge_factor = np.clip(charge_factor, 1 / self.charge_clip, self.charge_clip)
        
        return steric_cond * charge_factor
