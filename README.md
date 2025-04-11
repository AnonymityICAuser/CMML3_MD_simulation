# CMML3 Molecular Dynamics Simulation Repository

## Overview
This repository contains a complete workflow for molecular dynamics (MD) simulations of CMML3 systems, including analysis tools and visualization pipelines. The tools are optimized for GROMACS output files but can be adapted for other MD packages.


## Workflow Documentation

### 1. Simulation Setup
- Use GROMACS (or preferred MD package) with these recommended parameters:
  ```bash
  # Energy minimization
  emtol = 1000.0  # kJ/mol/nm
  emstep = 0.01   # nm

  # Equilibration
  nsteps = 50000   # 100ps
  dt = 0.002       # ps
  ```

### 2. GROMACS Analysis


### 3. Analysis Pipeline
1. Run `analysis_results.py` for comprehensive statistics
2. Generate specific comparisons with `Different_temps.py`
3. Create publication figures with `plot_xvg.py`



## Details in Analysis_visualization section

- Dependencies:
```text
Python >= 3.7
numpy >= 1.19
matplotlib >= 3.3
pandas >= 1.2
seaborn >= 0.11
scipy >= 1.6
```

#### 1. Different_temps.py
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
- Support for multi-column XVG files (e.g., per-residue RMSF)
- Configurable equilibration point detection (default: 50% of trajectory)

**Usage**:
```bash
# Basic usage
python analysis_results.py -i ./chains -o ./analysis_output

# Focused analysis (single temperature)
python analysis_results.py -i ./chains -o ./analysis_output -t 300K

# Focused analysis (single time point)
python analysis_results.py -i ./chains -o ./analysis_output -s 50ns
```

#### 3. plot_xvg.py
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

### chains/
**Directory Structure Convention**:
```
chains/
└── [ChainID]/               # e.g., ChainH, ChainD
    └── [Temperature]/       # e.g., 280K, 300K, 320K
        ├── [Timepoint]/     # e.g., 10ns, 50ns
        │   ├── rmsd.xvg
        │   ├── gyrate.xvg
        │   └── ...
        ├── temperature.xvg
        ├── pressure.xvg
        └── density.xvg
```

**Expected Analysis Files**:
1. Structural:
   - `rmsd.xvg`: Backbone RMSD
   - `rmsd_xtal.xvg`: Crystal reference RMSD
   - `rmsf.xvg`: Per-residue fluctuations
   - `gyrate.xvg`: Radius of gyration

2. Interactions:
   - `hb_all.xvg`: All hydrogen bonds
   - `hb_bb.xvg`: Backbone hydrogen bonds

3. Dynamics:
   - `eigenval.xvg`: PCA eigenvalues
   - `ev_components.xvg`: Eigenvector components
   - `ev_rmsf.xvg`: Eigenvector RMSF

### pyRMMA/
(For complete documentation, see [pyRMMA enhanced](https://github.com/yourusername/pyRMMA))

**Key Parameters**:
```python
# Typical configuration
parameters = {
    'correlation_times': [1e-9, 1e-8, 1e-7],  # Time scales for analysis
    'spectral_density': 'model-free',          # Options: 'model-free', 'lorentzian'
    'cutoff_frequency': 50.0,                  # MHz
    'relaxation_matrix': 'redfield'            # Matrix calculation method
}
```



