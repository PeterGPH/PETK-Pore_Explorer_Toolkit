# PETK GUI - Pore Explorer Toolkit for Nanopore Analysis

## Overview

PETK GUI (Pore Explorer Toolkit Graphical User Interface) is a comprehensive tool for analyzing molecular structures and their interactions with nanopores. The application provides an intuitive interface for molecular analysis, nanopore geometry design, and steric exclusion model (SEM) calculations.

## Features

### 🧬 Molecular Analysis (Tab 1: Nanopore Analyte Analysis)
- **PDB File Analysis**: Load and analyze protein structures from PDB files
- **Molecular Centering**: Automatically center molecules at the origin with surface alignment
- **Dimensional Analysis**: Calculate molecular dimensions, volumes, and bounding radii
- **Nanopore Compatibility**: Check if molecules fit through specified nanopore geometries
- **Quality Verification**: Verify centering and alignment quality with detailed scoring

### 🔬 SEM Calculation Setup (Tab 2: SEM Calculation Setup)
- **Pore Geometry Design**: Support for cylindrical and double-cone nanopore geometries
- **Parameter Configuration**: Set up simulation parameters including voltages, conductivities, and grid resolution
- **Python Environment Management**: Integrated conda environment setup for FEniCS calculations
- **Preview Generation**: Generate visualization previews of molecular movement through nanopores
- **Full SEM Simulation**: Run complete SEM calculations with current blockage analysis

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

### Tab 1: Molecular Analysis

#### Basic Workflow
1. **Set Project Information**
   - Select working directory
   - Configure nanopore parameters (diameter, thickness) for checking the fitness of analyte

2. **Load Molecular Structure**
   - Browse and select PDB file
   - Specify atom selection (default: "all")
   - Click "Analyze & Center" to process

3. **Review Results**
   - Check molecular dimensions and nanopore fit status
   - Review centering quality scores
   - Export analysis report if needed

4. **Visualization**
   - Use visualization controls to display molecule
   - Cycle through different representations
   - Center view as needed

#### Key Features
- **Automatic Centering**: Centers molecules at origin with optimal surface alignment
- **Fit Analysis**: Determines if molecules can pass through specified nanopores
- **Quality Scoring**: Provides detailed verification of centering accuracy
- **Export Options**: Save centered PDB files and analysis reports

### Tab 2: SEM Calculation Setup

#### Basic Workflow
1. **Select Pore Type**
   - Choose between cylindrical or double-cone geometries
   - Configure pore-specific parameters

2. **Configure Input Files**
   - Load centered analyte PDB (can sync from Tab 1)

3. **Set Simulation Parameters**
   - Define applied voltage and conductivities
   - Configure grid resolution and box dimensions
   - Set analyte movement parameters

4. **Environment Setup**
   - Create or select conda environment
   - Test Python/FEniCS installation
   - Validate all parameters

5. **Run Calculations**
   - Generate preview frames for visualization
   - Run full SEM simulation
  
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

### Output Files
- **Centered PDB**: Molecule centered at origin with surface alignment
- **JSON Config**: Configuration files for SEM calculation
- **Analysis Reports**: Detailed text reports of molecular analysis
- **SEM Results**: Current vs. position data from SEM calculations
- **Preview Images**: Visualization of molecular movement through nanopores

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

### Current Version 0.0.0
- Initial release with molecular analysis and SEM calculation capabilities assuming vertical passage of analyte through solid-state nanopore
- Support for cylindrical and double-cone pore geometries
- Integrated Python environment management
- Comprehensive validation and error checking

## Contact
If you have any feedback or encounter any bug, please contact Pinhao (Peter) via Pinhao2@illinois.edu