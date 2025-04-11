import os
import numpy as np
import matplotlib.pyplot as plt
import glob

# Base directory
base_dir = 'chains/ChainH'

# Temperature directories
temp_dirs = ['280K', '300K', '320K']
# Colors for each temperature
colors = {'280K': 'blue', '300K': 'green', '320K': 'red'}
# Time point directories
time_dirs = ['50ns'] # Using only 50ns data for better comparison

# Function to parse XVG files
def parse_xvg(filename):
    x = []
    y = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            data = line.strip().split()
            if len(data) >= 2:
                x.append(float(data[0]))
                # For gyrate.xvg, we only want the Rg value (column 1)
                if 'gyrate.xvg' in filename and len(data) > 1:
                    y.append(float(data[1]))  # Rg is in column 1
                else:
                    y.append(float(data[1]))  # Default to column 1 for other files
    return np.array(x), np.array(y)

# Create output directory for plots
output_dir = os.path.join(base_dir, 'joint_plots')
os.makedirs(output_dir, exist_ok=True)

# Dictionary to store file types
file_types = {
    'rmsd': 'rmsd.xvg',
    'rmsd_xtal': 'rmsd_xtal.xvg',
    'rmsf': 'rmsf.xvg',
    'gyrate': 'gyrate.xvg',
    'hb_all': 'hb_all.xvg',
    'hb_bb': 'hb_bb.xvg',
    'eigenval': 'eigenval.xvg',
    'ev_components': 'ev_components.xvg',
    'ev_rmsf': 'ev_rmsf.xvg',
    'num': 'num.xvg'
}

# Plot each file type separately
for file_type, file_name in file_types.items():
    plt.figure(figsize=(10, 6))
    
    for temp in temp_dirs:
        for time_dir in time_dirs:
            file_path = os.path.join(base_dir, temp, time_dir, file_name)
            
            if os.path.exists(file_path):
                x, y = parse_xvg(file_path)
                
                # Plot with appropriate label
                plt.plot(x, y, label=f'{temp}', color=colors[temp], linewidth=2)
    
    # Set plot title and labels
    plt.title(f'{file_type.upper()} at Different Temperatures (50ns)')
    plt.xlabel('Time (ps)')
    
    # Set y-label based on file type
    if file_type == 'rmsd' or file_type == 'rmsd_xtal':
        plt.ylabel('RMSD (nm)')
    elif file_type == 'rmsf':
        plt.ylabel('RMSF (nm)')
    elif file_type == 'gyrate':
        plt.ylabel('Radius of Gyration (nm)')
    elif file_type == 'hb_all' or file_type == 'hb_bb':
        plt.ylabel('Number of Hydrogen Bonds')
    elif file_type == 'eigenval':
        plt.ylabel('Eigenvalue')
    elif file_type == 'num':
        plt.ylabel('Number')
    else:
        plt.ylabel('Value')
    
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{file_type}_comparison.png'), dpi=300)
    plt.close()

# Also plot the temperature, pressure, and density files which are at the temperature directory level
for file_type in ['temperature', 'pressure', 'density']:
    plt.figure(figsize=(10, 6))
    
    for temp in temp_dirs:
        file_path = os.path.join(base_dir, temp, f'{file_type}.xvg')
        
        if os.path.exists(file_path):
            x, y = parse_xvg(file_path)
            
            # Plot with appropriate label
            plt.plot(x, y, label=f'{temp}', color=colors[temp], linewidth=2)
    
    # Set plot title and labels
    plt.title(f'{file_type.capitalize()} at Different Temperatures')
    plt.xlabel('Time (ps)')
    
    if file_type == 'temperature':
        plt.ylabel('Temperature (K)')
    elif file_type == 'pressure':
        plt.ylabel('Pressure (bar)')
    elif file_type == 'density':
        plt.ylabel('Density (kg/m³)')
    
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{file_type}_comparison.png'), dpi=300)
    plt.close()

print(f"All plots have been saved to {output_dir}")
