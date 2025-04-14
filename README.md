# CMML3 Molecular Dynamics Simulation: Scripts

## Information about my complex

| PDB ID | 2WWX |
|---------|-------|
| **Structure** | SidM/DrrA(GEF/GDF domain)-Rab1(GTPase) complex |
| **Resolution** | 1.50 Å |
| **Method** | X-RAY DIFFRACTION |
| **Organisms** | H. sapiens (Rab1), L. pneumophila (SidM/DrrA) |
| **Function** | Bacterial effector protein that hijacks host vesicular trafficking by acting as both a guanine nucleotide exchange factor (GEF) and GDI displacement factor (GDF) for Rab1 |
| **Significance** | Reveals mechanism by which L. pneumophila manipulates host membrane transport during infection |
| **Published** | EMBO J (2010) 29: 496, https://doi.org/10.1038/emboj.2009.347

## Overview
This repository contains a complete workflow for molecular dynamics (MD) simulations of CMML3 systems, including MD simulation step, analysis steps, and visualization pipelines. 

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
bash ./MD_sim_GROMACS/md_sim.bash INPUT_DIR OUTPUT_DIR
```
The final results in current work could be obtained from https://github.com/AnonymityICAuser/CMML3_additional_data. 

If you want to know more about the details, like actual parameters and setting in the code, please check the last section "Details in GROMACS simulation". 

### 2. GROMACS Analysis

Please check the ./MD_sim_GROMACS/analysis.bash. After the first simulation step, you could run the analysis directly, but please setting the correct working and output directories:

```bash
bash ./MD_sim_GROMACS/analysis.bash.
```

The result will be collected in folder `result_analysis/chains` for each chain and 'result_analysis/complex for the whole complex.

To reproduce the visualized results of ICA later, you could clone the repository of GROMACS analysis results:

```bash
git clone  https://github.com/AnonymityICAuser/CMML3_additional_data
```


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

##### 1. plot_xvg.py
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

##### 2. Different_temps.py
**Purpose**: Comparative analysis of simulation trajectories across multiple temperatures.

**Key Features**:
- Processes data from 280K, 300K, and 320K simulations (configurable)
- Generates 11 types of analysis plots:
  - Structural metrics: RMSD, RMSD (crystal reference)  - Geometric properties: Radius of gyration
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

Then run it directly:
```bash
python Different_temps.py
```


#### 3. analysis_results.py (optional)
**Purpose**: Comprehensive statistical comparative analysis and visualization of MD trajectories. The results didn't presented in final report.

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




##### ./pyRMMA
Python package for Ramachandran plot generation of protein structure, forked from https://github.com/gerdos/PyRAMA and has been modified. Please check details in https://github.com/AnonymityICAuser/PyRAMA_enhanced.

## Additional information: Details in GROMACS simulation 

### Structure Preparation
```bash
gmx pdb2gmx -f input.pdb -o complex_processed.gro -water spce -ignh -missing
```
Converts PDB to GROMACS format with SPCE water model.

### Box Definition
```bash
gmx editconf -f complex_processed.gro -o complex_box.gro -c -d 1.0 -bt cubic
```
Creates cubic simulation box with 1.0 nm protein-to-edge distance.

### Solvation
```bash
gmx solvate -cp complex_box.gro -cs spc216.gro -o complex_solv.gro -p topol.top
```
Fills box with water molecules.

### Ion Addition Setup
```bash
gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr -maxwarn 2
```
Prepares system for ion addition.

### Neutralization
```bash
gmx genion -s ions.tpr -o complex_ions.gro -p topol.top -neutral
```
Adds ions to neutralize system charge.

### Energy Minimization
```bash
gmx grompp -f em.mdp -c complex_ions.gro -p topol.top -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em -nb gpu -ntomp 8 -pin on
```
Removes unfavorable contacts and relaxes system.

### NVT Equilibration
```bash
gmx grompp -f NVT.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2
gmx mdrun -v -deffnm nvt -nb gpu -ntomp 8 -pin on
```
Constant Number, Volume, Temperature equilibration (200ps).

### NPT Equilibration
```bash
gmx grompp -f NPT.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2
gmx mdrun -v -deffnm npt -nb gpu -ntomp 8 -pin on
```
Constant Number, Pressure, Temperature equilibration (500ps).

### Production MD (1ns example)
```bash
gmx grompp -f md_1ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md_1ns.tpr -maxwarn 2
gmx mdrun -v -deffnm md_1ns -nb gpu -ntomp 8 -pin on
```
Production simulation run.

### Trajectory Processing
```bash
gmx trjconv -s md_1ns.tpr -f md_1ns.xtc -o md_1ns_noPBC.xtc -pbc mol -center
```
Removes PBC effects and centers protein for analysis.

## Temperature Variants

The pipeline runs parallel simulations at multiple temperatures (280K, 300K, 320K) by using temperature-specific parameter files for each stage.

## Simulation Extensions

The pipeline supports extending simulations to longer timescales (10ns, 50ns) by continuing from previous simulation endpoints:
```bash
gmx grompp -f md_10ns.mdp -c md_1ns.gro -t md_1ns.cpt -p topol.top -o md_10ns.tpr
gmx mdrun -v -deffnm md_10ns -nb gpu -ntomp 10 -pin on
```