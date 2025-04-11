import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def parse_xvg(file_path):
    x_data = []
    y_data = []
    
    title = ""
    x_label = ""
    y_label = ""
    labels = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('@'):
                parts = line.split()
                if len(parts) > 1:
                    if 'title' in parts[1].lower():
                        title = ' '.join(parts[2:]).replace('"', '')
                    elif 'xaxis' in parts[1].lower() and 'label' in parts[1].lower():
                        x_label = ' '.join(parts[2:]).replace('"', '')
                    elif 'yaxis' in parts[1].lower() and 'label' in parts[1].lower():
                        y_label = ' '.join(parts[2:]).replace('"', '')
                    elif 's' in parts[1].lower() and len(parts) > 2:
                        labels.append(' '.join(parts[2:]).replace('"', ''))
            
            elif not line.startswith('#') and not line.startswith('@'):
                values = line.split()
                if values:
                    try:
                        row_data = [float(val) for val in values]
                        if len(row_data) > 0:
                            x_data.append(row_data[0])
                            y_data.append(row_data[1:])
                    except ValueError:
                        continue

    x_data = np.array(x_data)
    y_data = np.array(y_data)
    
    return {
        'x_data': x_data,
        'y_data': y_data,
        'title': title,
        'x_label': x_label,
        'y_label': y_label,
        'labels': labels
    }

def plot_xvg(data, output_path=None):
    plt.figure(figsize=(10, 6))
    
    x_data = data['x_data']
    y_data = data['y_data']
    

    num_lines = y_data.shape[1] if len(y_data.shape) > 1 else 1
    
    if num_lines == 1:
        plt.plot(x_data, y_data, label=data['labels'][0] if data['labels'] else 'Data')
    else:
        for i in range(num_lines):
            label = data['labels'][i] if i < len(data['labels']) else f'Data {i+1}'
            plt.plot(x_data, y_data[:, i], label=label)
    
    if data['title']:
        plt.title(data['title'])
    if data['x_label']:
        plt.xlabel(data['x_label'])
    if data['y_label']:
        plt.ylabel(data['y_label'])
    
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, format='pdf', bbox_inches='tight')
        print(f"The figure has been saved in: {output_path}")
    else:
        plt.show()
    
    plt.close()

def process_nested_directories(source_dir, output_dir, pattern="*.xvg"):
    source_dir = Path(source_dir)
    output_dir = Path(output_dir)
    
    os.makedirs(output_dir, exist_ok=True)
    
    total_files = 0
    processed_files = 0
    
    for xvg_file in source_dir.glob('**/' + pattern):
        total_files += 1
        
        rel_path = xvg_file.relative_to(source_dir)
        output_path = output_dir / rel_path.with_suffix('.pdf')
        os.makedirs(output_path.parent, exist_ok=True)
        
        print(f"Processed: {xvg_file}")
        try:
            data = parse_xvg(xvg_file)

            means = np.mean(data['y_data'], axis=0)
            stds = np.std(data['y_data'], axis=0)
        
            print(f"File: {xvg_file.name}")
            for i, (mean, std) in enumerate(zip(means, stds)):
                label = data['labels'][i] if i < len(data['labels']) else f'Data {i+1}'
                print(f"  {label}: mean = {mean:.4f}, sd = {std:.4f}")
            
            plot_xvg(data, output_path)
            processed_files += 1
            
        except Exception as e:
            print(f"Error in {xvg_file}, {e}")
    
def main():
    parser = argparse.ArgumentParser(description='XVG Visualization Tool')
    parser.add_argument('-i', '--input' ,help='Containing directory of xvg files', required=True)
    parser.add_argument('-o', '--output', required=True, help='Output directory for PDF files')
    parser.add_argument('-p', '--pattern', default='*.xvg', help='Pattern to match xvg files (default: *.xvg)')
    
    args = parser.parse_args()
    
    process_nested_directories(args.input, args.output, args.pattern)

if __name__ == "__main__":
    main()