"""
Utility functions for SEM calculations.
Contains helper functions converted to DOLFINx.
"""

import numpy as np
import logging
import dolfinx
import dolfinx.fem as fem
import dolfinx.mesh as dmesh

# condfrac/readbinGrid moved to sem.grid_io (numpy-only, no dolfinx
# dependency) so pore_geometry/grid_io consumers can use them without
# pulling in dolfinx. Re-exported here for existing callers of sem.utils.
from .grid_io import condfrac, readbinGrid  # noqa: F401

logger = logging.getLogger(__name__)

def _safe_attr(obj, name):
    if obj is None or not hasattr(obj, name):
        return None
    try:
        value = getattr(obj, name)
        return value() if callable(value) else value
    except Exception:
        return None


def _is_dg0_space(V):
    element = V.ufl_element()

    family = _safe_attr(element, "family") or _safe_attr(element, "family_name")
    degree = _safe_attr(element, "degree")
    is_discontinuous = _safe_attr(element, "is_discontinuous") or _safe_attr(element, "discontinuous")

    basix_element = _safe_attr(element, "basix_element") or _safe_attr(element, "_basix_element")
    if basix_element is not None:
        if family is None:
            family = _safe_attr(basix_element, "family_name") or _safe_attr(basix_element, "family")
        if degree is None:
            degree = _safe_attr(basix_element, "degree")
        if is_discontinuous is None:
            is_discontinuous = _safe_attr(basix_element, "discontinuous")

    if degree != 0:
        return False

    if is_discontinuous is True:
        return True

    if family is None:
        return False

    family_str = str(family).lower()
    return "discontinuous" in family_str or family_str in ("dg", "dgp0", "dp")


def get_dof_coordinates(mesh_obj, V):
    """
    Return dof coordinates for a function space.
    Handles DG0 by using cell midpoints and the dofmap ordering.
    """
    if not _is_dg0_space(V):
        return V.tabulate_dof_coordinates()

    tdim = mesh_obj.topology.dim
    cell_map = mesh_obj.topology.index_map(tdim)
    num_cells = cell_map.size_local + cell_map.num_ghosts
    cells = np.arange(num_cells, dtype=np.int32)
    cell_midpoints = dmesh.compute_midpoints(mesh_obj, tdim, cells)

    dofmap = V.dofmap
    dof_index_map = dofmap.index_map
    num_dofs = dof_index_map.size_local + dof_index_map.num_ghosts
    coords = np.zeros((num_dofs, mesh_obj.geometry.dim), dtype=cell_midpoints.dtype)

    cell_dofs = dofmap.list
    if hasattr(cell_dofs, "links"):
        for cell in range(num_cells):
            dofs = cell_dofs.links(cell)
            if len(dofs) != 1:
                raise ValueError("DG0 space is expected to have exactly one dof per cell")
            coords[dofs[0]] = cell_midpoints[cell]
    else:
        cell_dofs = np.asarray(cell_dofs)
        if cell_dofs.ndim != 2 or cell_dofs.shape[1] != 1:
            raise ValueError("DG0 space is expected to have exactly one dof per cell")
        coords[cell_dofs[:, 0]] = cell_midpoints

    return coords


def loadFunc(mesh_obj, V, sig_func, interpfunction, bulk_conductivity):
    """
    Load function values using interpolation (DOLFINx version).
    Enhanced with better error handling and debugging.
    """
    try:
        logger.info("Starting loadFunc (DOLFINx version)...")
        
        # Get DOF coordinates (DOLFINx way)
        x = get_dof_coordinates(mesh_obj, V)
        logger.info(f"DOF coordinates shape: {x.shape}")
        
        # Evaluate interpolation function at DOF coordinates
        values = interpfunction(x)
        logger.info(f"Interpolated {len(values)} values")
        
        # Handle NaN values if any
        nan_mask = np.isnan(values)
        if np.any(nan_mask):
            logger.warning(f"Found {np.sum(nan_mask)} NaN values, replacing with bulk conductivity")
            # Assume bulk conductivity is the fill_value from the interpolator
            # or use a reasonable default
            values[nan_mask] = bulk_conductivity
        
        # Set values in function (DOLFINx way)
        sig_func.x.array[:] = values
        sig_func.x.scatter_forward()
        
        logger.info("loadFunc completed successfully (DOLFINx)")
        
    except Exception as e:
        logger.error(f"Error in loadFunc (DOLFINx): {e}")
        logger.error(f"Error type: {type(e)}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        raise
