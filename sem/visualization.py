"""
Visualization functions for SEM calculations.
Provides plotting and preview capabilities.
"""

import numpy as np
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def plot_conductivity_slice(x_coords, z_coords, conductivity_2d, z_position, 
                                 pore_type, pore_params=None, save_path=None, show=True):
    """
    Plot conductivity using pcolormesh - much better coordinate control than imshow.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    print(f"🎨 Using pcolormesh for better coordinate control")
    
    # Create coordinate meshgrids
    X_mesh, Z_mesh = np.meshgrid(x_coords, z_coords, indexing='ij')
    
    print(f"Mesh shapes: Z_mesh {Z_mesh.shape}, X_mesh {X_mesh.shape}, data {conductivity_2d.shape}")
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Use pcolormesh - much more accurate coordinate mapping
    pc = plt.pcolormesh(X_mesh, Z_mesh, conductivity_2d,
                        cmap="viridis", shading="auto")

    cb = plt.colorbar(pc)
    cb.set_label("Conductivity (S/m)", fontsize=16, fontweight="bold")
    cb.ax.tick_params(labelsize=12)  # colorbar tick font size

    plt.xlabel("X position (Å)", fontsize=14, fontweight="bold")
    plt.ylabel("Z position (Å)", fontsize=14, fontweight="bold")
    # Add membrane boundaries
    if pore_params:
        membrane_thickness = pore_params.get('membrane_thickness', 200)
        membrane_z_offset = pore_params.get('membrane_z_offset', 0) or 0.0
        membrane_half = membrane_thickness / 2
        
        # These will be EXACTLY at the right coordinates
        plt.axhline(y=membrane_half+membrane_z_offset, color='white', linestyle='-', linewidth=3, alpha=0.9,
                   label=f'Membrane boundaries (±{membrane_half:.1f}Å)')
        plt.axhline(y=-membrane_half+membrane_z_offset, color='white', linestyle='-', linewidth=3, alpha=0.9)
        
        # Add membrane region highlighting
        plt.axhspan(-membrane_half+membrane_z_offset, membrane_half+membrane_z_offset, color='red', alpha=0.1, 
                   label=f'Membrane region ({membrane_thickness}Å thick)')
    
    # Add analyte position
    plt.axhline(y=z_position, color='yellow', linestyle='--', linewidth=2, alpha=0.8,
               label=f'Analyte at Z = {z_position:.1f}Å')
    
    plt.title(f'Conductivity Map - {pore_type.title()} Pore (pcolormesh)')
    # plt.legend()
    plt.grid(True, alpha=0.3)
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Saved pcolormesh plot: {save_path}")
    
    if show:
        plt.show()
    
    return True

def plot_results(results, save_path=None, show=True):
    """
    Plot SEM simulation results.
    
    Args:
        results: Results dictionary from SEM simulation
        save_path: Path to save plot (optional)
        show: Whether to show the plot
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("Matplotlib not available, cannot create plots")
        return None
    
    z_positions = results['z_positions']
    currents = results['currents']
    normalized_currents = results['normalized_currents']
    blockages = results['blockages']
    
    # Create subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    
    # Current vs position
    ax1.plot(z_positions, currents * 1e12, 'b-', linewidth=2)
    ax1.set_xlabel('Z Position (Å)')
    ax1.set_ylabel('Current (pA)')
    ax1.set_title(f'Current vs Z Position - {results["pore_type"].title()} Pore')
    ax1.grid(True, alpha=0.3)
    
    # Normalized current vs position
    ax2.plot(z_positions, normalized_currents, 'g-', linewidth=2)
    ax2.set_xlabel('Z Position (Å)')
    ax2.set_ylabel('Normalized Current')
    ax2.set_title('Normalized Current vs Z Position')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1.1)
    
    # Blockage vs position
    ax3.plot(z_positions, blockages, 'r-', linewidth=2)
    ax3.set_xlabel('Z Position (Å)')
    ax3.set_ylabel('Blockage (%)')
    ax3.set_title('Current Blockage vs Z Position')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        logger.info(f"Saved results plot: {save_path}")
    
    if show:
        plt.show()
    
    plt.close()
    
    return True

def create_preview_frames(sem_instance, num_frames=4, save_plots=True, output_dir="."):
    """
    Create preview frames showing analyte movement.
    
    Args:
        sem_instance: VerticalMovementSEM instance
        num_frames: Number of preview frames to create
        save_plots: Whether to save plots
        output_dir: Directory to save plots
    """
    logger.info(f"Creating {num_frames} preview frames of analyte movement...")
    logger.info(f"Pore type: {sem_instance.pore_type}")
    
    # Calculate total number of steps
    total_steps = int(abs(sem_instance.z_end - sem_instance.z_start) / sem_instance.z_step) + 1
    z_positions = np.linspace(sem_instance.z_start, sem_instance.z_end, total_steps)
    
    logger.info(f"Total available steps: {total_steps}")
    logger.info(f"Selecting {num_frames} frames from {total_steps} total steps")
    
    # Select frames to preview
    if num_frames >= total_steps:
        frame_indices = list(range(total_steps))
        logger.info(f"Using all {total_steps} available positions")
    else:
        frame_indices = np.linspace(0, total_steps-1, num_frames, dtype=int)
        selected_z_positions = [z_positions[idx] for idx in frame_indices]
        logger.info(f"Selected Z positions: {[f'{z:.1f}' for z in selected_z_positions]}")
    
    # Collect pore parameters for plotting
    pore_params = {
        'pore_radius': sem_instance.pore_radius,
        'membrane_thickness': sem_instance.membrane_thickness,
        'membrane_z_offset': sem_instance.membrane_z_offset
    }
    
    if sem_instance.pore_type == "cylindrical":
        pore_params['corner_radius'] = sem_instance.corner_radius
    elif sem_instance.pore_type == "double_cone":
        pore_params['outer_radius'] = sem_instance.outer_radius
    elif sem_instance.pore_type == "biological":
        pore_params['biological_pore_pdb'] = sem_instance.biological_pore_pdb
    elif sem_instance.pore_type == "bin_file":
        pore_params['bin_file_path'] = sem_instance.bin_file_path
    
    for i, frame in enumerate(frame_indices):
        z_pos = z_positions[frame]
        logger.info(f"Creating preview frame {i+1}/{len(frame_indices)} at Z = {z_pos:.1f} Å")
        
        # Get conductivity grid for this position
        x_coords, z_coords, cond = sem_instance.get_conductivity_grid_for_preview(z_pos)
        
        # Create save path if needed
        save_path = None
        if save_plots:
            save_path = Path(output_dir) / f'{sem_instance.output_prefix}_{sem_instance.pore_type}_preview_frame_{frame:04d}_z_{z_pos:.1f}A.png'
        
        # Plot the frame
        plot_conductivity_slice(
            x_coords, z_coords, cond, z_pos,
            sem_instance.pore_type, pore_params,
            save_path=save_path, show=not save_plots
        )
        
        # Also save data if requested
        if save_plots:
            data_path = Path(output_dir) / f'{sem_instance.output_prefix}_{sem_instance.pore_type}_preview_frame_{frame:04d}_z_{z_pos:.1f}A.dat'
            np.savetxt(data_path, cond)
    
    logger.info("Preview frames created successfully!")

def export_dx_file(sem_instance, filename=None, custom_resolution=None):
    """
    Simple, direct export to DX format - mimics the gridData workflow pattern.
    Just reads the data and exports it without complex parameter handling.
    """
    try:
        from gridData import Grid
    except ImportError:
        logger.error("gridData library not found. Please install with: pip install gridData")
        return None

    if filename is None:
        filename = f"{sem_instance.output_prefix}_{sem_instance.pore_type}_pore_conductivity.dx"

    logger.info(f"Writing pore conductivity map to DX format: {filename}")

    # Unified sampling using mesh domain attributes
    if hasattr(sem_instance, 'domain_min') and hasattr(sem_instance, 'domain_max') and hasattr(sem_instance, 'num_cells'):
        logger.info("🚀 Sampling interpolator using DOLFINx mesh domain...")

        if custom_resolution is not None:
            # Compute number of cells based on custom resolution
            num_cells = [
                int((sem_instance.domain_max[i] - sem_instance.domain_min[i]) / custom_resolution)
                for i in range(3)
            ]
            delta = [custom_resolution] * 3
            logger.info(f"Using custom resolution: {custom_resolution} Å")
            logger.info(f"Computed num_cells: {num_cells}")
        else:
            num_cells = sem_instance.num_cells
            delta = [
                (sem_instance.domain_max[i] - sem_instance.domain_min[i]) / num_cells[i]
                for i in range(3)
            ]
            logger.info(f"Using mesh resolution: {delta}")

        num_points = [n + 1 for n in num_cells]

        x_range = np.linspace(sem_instance.domain_min[0], sem_instance.domain_max[0], num_points[0])
        y_range = np.linspace(sem_instance.domain_min[1], sem_instance.domain_max[1], num_points[1])
        z_range = np.linspace(sem_instance.domain_min[2], sem_instance.domain_max[2], num_points[2])

        nx, ny, nz = len(x_range), len(y_range), len(z_range)

        logger.info(f"📊 Sampling interpolator: {nx} × {ny} × {nz} grid")

        # Create meshgrid and sample
        X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing='ij')
        grid_coords = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

        # Sample the interpolator
        conductivity_flat = sem_instance.base_cond_interp(grid_coords)
        conductivity_data = conductivity_flat.reshape((nx, ny, nz))

        # Handle NaN values
        nan_mask = np.isnan(conductivity_data)
        if np.any(nan_mask):
            conductivity_data[nan_mask] = sem_instance.bulk_conductivity
            logger.info(f"⚠️  Replaced {np.sum(nan_mask):,} NaN values")

        origin = sem_instance.domain_min

        logger.info(f"✅ Interpolator sampled successfully")
        logger.info(f"   Data range: [{np.min(conductivity_data):.6f}, {np.max(conductivity_data):.6f}] S/m")
        logger.info(f"   Delta: {delta}")
        logger.info(f"   Origin: {origin}")

    else:
        logger.error("❌ Mesh domain attributes not available for export")
        return None

    # Export
    logger.info("💾 Creating Grid object...")
    
    try:
        pore_grid = Grid(
            grid=conductivity_data,
            delta=delta,
            origin=origin
        )

        pore_grid.export(filename)

        logger.info(f"✅ DX file written successfully: {filename}")
        logger.info(f"📊 Final export info:")
        logger.info(f"   File: {filename}")
        logger.info(f"   Shape: {conductivity_data.shape}")
        logger.info(f"   Delta: {delta}")
        logger.info(f"   Origin: {origin}")
        logger.info(f"   Data range: [{np.min(conductivity_data):.6f}, {np.max(conductivity_data):.6f}] S/m")
        logger.info(f"🚀 VMD command: vmd {filename}")

        return filename
        
    except Exception as e:
        logger.error(f"❌ Failed to create or export DX file: {e}")
        return None


def export_mesh(sem_instance, filename=None):
    """
    Export the DOLFINx mesh to XDMF for visualization (e.g., ParaView).
    """
    try:
        import dolfinx.io as io
    except ImportError:
        logger.error("dolfinx is not available; cannot export mesh.")
        return None

    if filename is None:
        filename = f"{sem_instance.output_prefix}_{sem_instance.pore_type}_mesh.xdmf"

    logger.info("Writing mesh to XDMF: %s", filename)
    with io.XDMFFile(sem_instance.comm, filename, "w") as xdmf:
        xdmf.write_mesh(sem_instance.mesh)

    return filename
