# CMML3 Molecular Dynamics Simulation Repository

## Overview
This repository contains a complete workflow for molecular dynamics (MD) simulations of CMML3 systems, including analysis steps and visualization pipelines. The tools are optimized for GROMACS output files but can be adapted for other MD packages.


## Workflow Documentation

### 1. GROMACS MD Simulation 
Please check the ./MD_sim_GROMACS/md_sim.bash and relevant setting files .mdp. 

Here presents the meaning of each type of mdp file:
- em.mdp: Energy minimization parameter file. Used to relax the initial structure and remove steric clashes before running dynamics.
- ions.mdp: Parameter file for ion placement. Used when adding ions to neutralize the system or achieve a specific ionic concentration.
- NVT_mdps folder: Contains parameter files for NVT ensemble simulations (constant Number of particles, Volume, and Temperature) at different temperatures. NVT equilibration is typically the first equilibration phase after energy minimization.
- NPT_mdps folder: Contains parameter files for NPT ensemble simulations (constant Number of particles, Pressure, and Temperature) at different temperatures. NPT equilibration typically follows NVT equilibration to allow the system volume to adjust.
- Production MD files:
	1. md_50ns_280K.mdp: 50 nanosecond production run at 280K
	2. md_50ns_300K.mdp: 50 nanosecond production run at 300K
	3. md_50ns_320K.mdp: 50 nanosecond production run at 320K
These are for the actual data collection simulations after equilibration.

Run the bash script directly when you are ready for simulation:
```bash 
bash ./MD_sim_GROMACS/md_sim.bash
```


### 2. GROMACS Analysis
Please check the ./MD_sim_GROMACS/analysis.bash. After the first simulation step, you could run the analysis directly, but please setting the correct working and output direcotries:

```bash
bash ./MD_sim_GROMACS/analysis.bash.
```

The result will be collected in folder `result_analysis/chains` for each chain and 'result_analysis/complex for the whole complex.


### 3. Visualization Pipeline

The visualization is based on Python.

- Dependencies:
```text
Python >= 3.7
numpy >= 1.19
matplotlib >= 3.3
pandas >= 1.2
seaborn >= 0.11
scipy >= 1.6
```

##### 1. Different_temps.py
**Purpose**: Comparative analysis of simulation trajectories across multiple temperatures.

**Key Features**:
- Processes data from 280K, 300K, and 320K simulations (configurable)
- Generates 11 types of analysis plots:
  - Structural metrics: RMSD, RMSD (crystal reference), RMSF
  - Geometric properties: Radius of gyration
  - Interaction analysis: Hydrogen bonds (all and backbone-only)
  - Dynamic modes: Eigenvalues, eigenvector components, eigenvector RMSF
  - System properties: Temperature, pressure, density
- Outputs PNG images with standardized formatting (300 DPI)

**Configuration**:
```python
# Base directory containing simulation data
base_dir = 'chains/ChainH'  # Modify to match your data structure

# Temperature directories to analyze
temp_dirs = ['280K', '300K', '320K']  # Add/remove temperatures as needed

# Time points to analyze (currently only 50ns)
time_dirs = ['50ns']  # Can extend to ['10ns', '50ns', '100ns']
```

#### 2. analysis_results.py
**Purpose**: Comprehensive statistical analysis and visualization of MD trajectories.

**Analysis Capabilities**:
- Automatic recursive discovery of XVG files in directory trees
- Multi-dimensional analysis:
  - By temperature (e.g., 280K vs 300K vs 320K)
  - By simulation time (e.g., 10ns vs 50ns)
  - By metric type (RMSD, RMSF, etc.)
- Statistical outputs:
  - Full-trajectory and equilibrated-region means
  - Standard deviations
  - Comparison heatmaps

**Advanced Features**:
- Robust error handling for malformed XVG files
- Automatic label extraction from XVG metadata
- Support for multi-column XVG files 

**Usage**:
```bash
# Basic usage
python analysis_results.py -i ./chains -o ./analysis_output

# Focused analysis (single temperature)
python analysis_results.py -i ./chains -o ./analysis_output -t 300K

# Focused analysis (single time point)
python analysis_results.py -i ./chains -o ./analysis_output -s 50ns
```

##### 3. plot_xvg.py
**Purpose**: Quick visualization of individual XVG files with publication-quality output.

**Features**:
- Preserves original directory structure in output
- Automatic extraction of titles and axis labels from XVG headers
- Statistical summary (mean ± SD) printed to console
- Configurable file patterns (default: *.xvg)

**Usage**:
```bash
# Process all XVG files recursively
python plot_xvg.py -i ./raw_data -o ./plots

# Process specific file patterns
python plot_xvg.py -i ./raw_data -o ./plots -p "*rmsd*.xvg"
```

**Expected Analysis Files**:
1. Structural:
   - `rmsd.xvg`: Backbone RMSD
   - `rmsd_xtal.xvg`: Crystal reference RMSD
   - `gyrate.xvg`: Radius of gyration

2. Interactions:
   - `hb_all.xvg`: All hydrogen bonds
   - `hb_bb.xvg`: Backbone hydrogen bonds


#####./pyRMMA
Python package for Ramachandran plot generation of protein structure, forked from https://github.com/gerdos/PyRAMA and has been modified. Please check details in https://github.com/AnonymityICAuser/PyRAMA_enhanced.


