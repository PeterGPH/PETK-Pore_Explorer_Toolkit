# PETK GUI - Pore Explorer Toolkit for Nanopore Analysis

## Overview

PETK GUI (Pore Explorer Toolkit Graphical User Interface) is a comprehensive tool for analyzing molecular structures and their interactions with nanopores. The application provides an intuitive interface for molecular analysis, nanopore geometry design, and steric exclusion model (SEM) calculations.

## Installation

### Prerequisites
- **VMD**: For molecular visualization and analysis
- **Python 3.9+**: For SEM calculations
- **Conda**: For Python environment management

### Required Python Packages
The application can automatically create a conda environment with these packages, which can be installed easily via one click on button "Create SEM Env" on tab2 :
- FEniCS (2019.1.0)
- NumPy
- SciPy
- Matplotlib
- MDAnalysis
- h5py

### Setup Instructions

1. **Clone or download the PETK GUI files**
   ```bash
   # Ensure you have these files:
   # - petk_gui.tcl (main GUI script)
   # - vertical_movement_sem.py (Python SEM calculation module)
   ```

2. **Install VMD** (if not already installed)
   - Download from: https://www.ks.uiuc.edu/Research/vmd/

3. **Install Conda** (if not already installed)
   - Download from: https://conda.io/miniconda.html

4. **Launch the application**
   ```bash
   # From VMD Tk Console or command line with VMD in path:
   source petk_gui.tcl
   ::PETK::gui::petk_gui
   ```

## Usage Guide

### Tab 1: Nanopore Setup

#### Basic Workflow
1. **Configure Output Settings**
   - Select working directory
   - Set output file prefix

2. **Set Box Dimensions**
   - Choose auto-calculated (recommended) or manual mode
   - Configure distance cutoff for auto-calculation
   - View real-time dimension calculations

3. **Design Pore Geometry**
   - Select pore type: Cylindrical or Double-Cone
   - Configure pore-specific parameters:
     - **Cylindrical**: Diameter, corner radius, thickness
     - **Double-Cone**: Inner diameter, outer diameter, thickness

4. **Visualize and Validate**
   - View real-time 3D pore visualization in VMD
   - Check parameter validity status
   - Update preview as needed

#### Key Features
- **Native VMD Integration**: Real-time 3D visualization with membrane and pore geometry
- **Automatic Validation**: Comprehensive checking of pore-box compatibility
- **Visual Feedback**: Color-coded status indicators for invalid configurations

### Tab 2: Analyte Setup

#### Basic Workflow
1. **Load Molecular Structure**
   - Browse and select PDB file
   - Specify atom selection (default: "all")
   - Click "Analyze & Center" to process

2. **Review Analysis Results**
   - Check molecular dimensions and nanopore fit status
   - Review centering and alignment quality scores
   - View detailed molecular information

3. **Visualization Controls**
   - Show/hide molecule in VMD
   - Change molecular representations
   - Center view and cycle through display styles

4. **Export Results**
   - Export detailed analysis report
   - Save centered PDB file for SEM calculations

#### Enhanced Features
- **Intelligent Surface Alignment**: Automatically orients molecules for optimal pore passage
- **Quality Scoring**: 12-point verification system for centering accuracy
- **Fit Analysis**: Real-time compatibility checking with designed nanopore
- **Detailed Reporting**: Comprehensive molecular analysis with exportable reports

### Tab 3: SEM Setup

#### Basic Workflow
1. **Configure Simulation Parameters**
   - Set applied voltage and conductivities
   - Configure grid resolution and VdW radii settings
   - Define analyte movement parameters (Z range and step size)

2. **Python Environment Setup**
   - Create or select conda environment
   - Test Python/FEniCS installation using "Test Python" button
   - Use "Create SEM Env" for automatic environment setup

3. **Validate and Run**
   - Use "Validate Parameters" to check all settings
   - Generate trajectory preview with "Preview Simulation"
   - Run full SEM calculation with "Run Simulation"

#### Key Features
- **One-Click Environment Setup**: Automatic creation of complete FEniCS environment
- **Real-time Validation**: Comprehensive parameter checking before calculation
- **Trajectory Preview**: Visual preview of analyte movement through nanopore
- **Grid Estimation**: Real-time calculation of computational requirements

## Enhanced Features

### Native VMD Visualization
- **Real-time 3D Pore Rendering**: Live visualization of nanopore geometry
- **Membrane Representation**: Visual membrane with pore cutouts
- **Parameter Integration**: Automatic updates when parameters change
- **Validation Feedback**: Visual warnings for invalid configurations

  
## Python SEM Module

The included `vertical_movement_sem.py` module provides:
- **FEniCS-based SEM calculations**
- **Van der Waals radii database**
- **Multiple pore geometries** (cylindrical, double-cone)
- **JSON configuration support**
- **Visualization capabilities**
- **Current blockage analysis**

### Command Line Usage
```bash
# Run full simulation
python vertical_movement_sem.py config.json run

# Generate preview only
python vertical_movement_sem.py config.json preview_only
```

## File Formats

### Input Files
- **PDB Files**: Protein Data Bank format for molecular structures
- 
### Output Files
- **Centered PDB**: Molecule centered at origin with surface alignment
- **JSON Config**: Configuration files for SEM calculation
- **Analysis Reports**: Detailed text reports of molecular analysis
- **SEM Results**: Current vs. position data from SEM calculations
- **Trajectory Files**: DCD trajectories for movement preview
- **Preview Images**: Visualization of molecular movement through nanopores (optional)

## Troubleshooting

### Common Issues

**GUI Won't Launch**
- Ensure VMD is properly installed and in PATH
- Check that Tcl/Tk is available
- Verify script permissions

**Python Environment Issues**
- Use "Create SEM Env" button to set up environment automatically
- Ensure conda is in PATH
- Check FEniCS installation with "Test Python" button

**SEM Calculation Errors**
- Validate all parameters before running
- Check that analyte PDB file exists and is readable
- Ensure sufficient computational resources for grid resolution (3.4M - 3 min, 27M - 15 min)

**Memory Issues**
- Reduce grid resolution for large systems
- Use auto-calculated box dimensions
- Consider reducing number of simulation steps


## Configuration Examples

### Basic Cylindrical Pore
```json
{
  "pore_geometry": {
    "pore_type": "cylindrical",
    "pore_radius": 100.0,
    "corner_radius":50.0,
    "membrane_thickness": 200.0
  }
}
```

### Double-Cone Pore
```json
{
  "pore_geometry": {
    "pore_type": "double_cone",
    "pore_radius": 75.0,
    "outer_radius": 150.0,
    "membrane_thickness": 200.0
  }
}
```

## Technical Details

### Units
- All distances in Angstroms (Å)
- Voltages in millivolts (mV)
- Conductivities in Siemens per meter (S/m)
- Currents in Amperes (nA)

### Algorithms
- **Centering**: Uses center of mass with surface normal alignment
- **SEM Calculations**: Finite element method with FEniCS
- **Conductivity Modeling**: Distance-based conductivity reduction
- **Van der Waals Radii**: Element-specific atomic radii from literature

## Citation

If you use PETK GUI in your research, please cite the original paper:
```
@article{wilson2019rapid, title={Rapid and accurate determination of nanopore ionic current using a steric exclusion model},
author={Wilson, James and Sarthak, Kumar and Si, Wei and Gao, Luyu and Aksimentiev, Aleksei},
journal={Acs Sensors},
volume={4},
number={3},
pages={634--644},
year={2019},
publisher={ACS Publications}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Version History
### Current Version 0.0.1

- New 3-tab interface: Separated nanopore design, analyte setup, and SEM calculation
- Native VMD visualization: Real-time 3D pore geometry rendering
- Trajectory preview: Visual preview of analyte movement through nanopores



### Previous Version 0.0.0
- Initial release with molecular analysis and SEM calculation capabilities assuming vertical passage of analyte through solid-state nanopore
- Support for cylindrical and double-cone pore geometries
- Integrated Python environment management
- Comprehensive validation and error checking

## Contact
If you have any feedback or encounter any bug, please contact Pinhao (Peter) via Pinhao2@illinois.edu
