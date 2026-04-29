"""
Utility functions for SEM calculations.
Contains helper functions converted to DOLFINx.
"""

import os
import numpy as np
import logging
import dolfinx
import dolfinx.fem as fem
import dolfinx.mesh as dmesh

logger = logging.getLogger(__name__)

def condfrac(invec):
    """
    Convert values to conductivity fractions (from original code).
    Points on line (min,0) (max,1)
    """
    minr = 1.3
    maxr = 4.1
    slope = 1.0/(maxr-minr)
    int_val = -minr*slope
    result = slope*invec + int_val
    result[result<0] = 0.0000001
    result[result>1] = 1.0
    return result

def readbinGrid(name, mask_radius=-1, *, return_metadata=False):
    """
    Read binary grid file (from original code).
    
    Args:
        name: Path to binary grid file.
        mask_radius: Optional radius for masking values.
        return_metadata: If True, return a metadata dict with origin, spacing, and grid shape.
    
    Returns:
        Tuple containing the 3D values, physical dimensions, grid counts,
        and optionally metadata about the grid spacing/origin.
    """
    if not os.path.isfile(name):
        print(name+" doesn't exist, EXITING")
        exit()
    
    with open(name, 'rb') as f:
        val1d = np.fromfile(f, dtype=np.float32)

    if val1d.size < 7:
        raise ValueError(f"Binary grid file {name} is too small to contain a header")

    resolution = float(val1d[6])
    if resolution <= 0:
        raise ValueError(f"Invalid grid spacing ({resolution}) recorded in {name}")

    delta = np.array([resolution, resolution, resolution], dtype=np.float32)
    origin = np.array([val1d[3], val1d[4], val1d[5]], dtype=np.float32)
    shape = (
        int(np.ceil(val1d[0])),
        int(np.ceil(val1d[1])),
        int(np.ceil(val1d[2]))
    )

    expected_values = shape[0] * shape[1] * shape[2]
    data = val1d[7:]

    if data.size != expected_values:
        raise ValueError(
            f"Binary grid {name} contains {data.size} values but expected {expected_values}"
        )

    val3d = np.reshape(data, shape, order='F')
    
    if mask_radius>0:
        x_ = np.arange(origin[0],origin[0]+shape[0]*delta[0],delta[0])
        y_ = np.arange(origin[1],origin[1]+shape[1]*delta[1],delta[1])
        z_ = np.arange(origin[2],origin[2]+shape[2]*delta[2],delta[2])
        assert len(x_) == val3d.shape[0], "x is wrong size"
        assert len(y_) == val3d.shape[1], "y is wrong size"
        assert len(z_) == val3d.shape[2], "z is wrong size"
        xx,yy,zz = np.meshgrid(x_,y_,z_, indexing='ij')

        msk = xx*xx+yy*yy>mask_radius*mask_radius
        val3d[msk] = 0.00001
        
    L = delta[0]*shape[0]; W = delta[1]*shape[1]; H = delta[2]*shape[2]
    nx = int(shape[0])
    ny = int(shape[1])
    nz = int(shape[2])
    Lm = L-delta[0]
    Wm = W-delta[1]
    Hm = H-delta[2]

    if return_metadata:
        metadata = {
            "origin": origin,
            "spacing": delta,
            "grid_shape": shape,
            "resolution": resolution,
        }
        return val3d, [Lm,Wm,Hm], [nx,ny,nz], metadata

    return val3d, [Lm,Wm,Hm], [nx,ny,nz]

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
