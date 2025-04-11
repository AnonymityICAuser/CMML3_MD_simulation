import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import pandas as pd
from scipy import stats
import seaborn as sns
from collections import defaultdict
import traceback # For better error reporting

def parse_xvg(file_path):
    """Parse XVG file and return data and metadata"""
    x_data = []
    y_data = []

    title = ""
    x_label = ""
    y_label = ""
    labels = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('@'):
                    # Use maxsplit to handle potential spaces in labels/titles
                    parts = line.split(maxsplit=2)
                    if len(parts) > 1:
                        command = parts[1].lower()
                        # Ensure value exists before trying to replace quotes
                        value = parts[2].replace('"', '') if len(parts) > 2 else ""

                        if 'title' in command:
                            title = value
                        elif 'xaxis' in command and 'label' in command:
                            x_label = value
                        elif 'yaxis' in command and 'label' in command:
                            y_label = value
                        # Improved label parsing: Check for 's' followed by a number
                        elif command.startswith('s') and command[1:].isdigit() and value:
                             labels.append(value)

                elif not line.startswith('#') and not line.startswith('@'):
                    values = line.split()
                    if values:
                        try:
                            # Ensure all values can be converted to float
                            row_data = [float(val) for val in values]
                            if not row_data: # Skip empty lines after split
                                continue
                            x_data.append(row_data[0])
                            # Handle cases with only x column (though less common)
                            # Append an empty list if no y-values, will become np.nan later
                            y_data.append(row_data[1:] if len(row_data) > 1 else [])
                        except ValueError:
                            # Skip lines with non-numeric data that aren't comments/headers
                            # print(f"Skipping non-numeric data line in {file_path}: {line}")
                            continue

        x_data = np.array(x_data, dtype=float)

        # Handle potentially jagged y_data (different numbers of columns per row)
        # Find the maximum number of y-columns
        max_cols = 0
        if y_data:
            max_cols = max(len(row) for row in y_data)

        # Pad shorter rows with NaN
        y_data_padded = []
        for row in y_data:
            padded_row = row + [np.nan] * (max_cols - len(row))
            y_data_padded.append(padded_row)

        # Convert to numpy array, handling the case where y_data might be empty
        y_data = np.array(y_data_padded, dtype=float) if y_data_padded else np.array([])

        # Ensure y_data is 2D, even if only one column or empty
        if y_data.ndim == 1 and len(x_data) > 0:
             y_data = y_data.reshape(-1, 1)
        elif len(x_data) == 0: # If no data rows were found
             y_data = np.empty((0, max(1, len(labels)))) # Use labels count or 1 col

    except FileNotFoundError:
        print(f"Error: File not found {file_path}")
        return None
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        traceback.print_exc() # Print detailed traceback
        return None

    return {
        'x_data': x_data,
        'y_data': y_data,
        'title': title,
        'x_label': x_label,
        'y_label': y_label,
        'labels': labels if labels else [y_label] if y_label else ['Y-Value'] # Provide default labels
    }

def collect_all_data(root_dir):
    """
    Collect data from all XVG files recursively, grouping by path components.
    Assumes structure like: ... / group_id / tempK / time_ns / analysis_type.xvg
    """
    root_dir = Path(root_dir).resolve() # Ensure absolute path
    # Use defaultdict for easier nested dictionary creation
    all_data = defaultdict(lambda: defaultdict(dict))

    print(f"Searching for XVG files in: {root_dir}")
    found_files = list(root_dir.rglob('*.xvg'))
    print(f"Found {len(found_files)} potential XVG files.")

    if not found_files:
         print("Warning: No .xvg files found recursively.")
         return {} # Return empty dict if no files found

    for xvg_file in found_files:
        try:
            # Extract relevant parts from the path
            parts = xvg_file.parts
            # Find the index of the root directory to slice relative path
            try:
                root_index = parts.index(root_dir.name)
                relative_parts = parts[root_index + 1:]
            except ValueError:
                print(f"Warning: Could not determine relative path for {xvg_file}. Skipping.")
                continue

            # Flexible path parsing - adjust indices based on your *actual* structure depth
            # Example: root / category / group / temp / time / file.xvg
            #          idx 0 / idx 1 / idx 2 / idx 3 / idx 4 / idx 5 (relative to root)
            # We expect temp and time directories just before the file.
            # The directory before temp is assumed to be the group_id.
            if len(relative_parts) >= 4: # Need at least group/temp/time/file.xvg
                time_str = relative_parts[-2]
                temp_str = relative_parts[-3]
                group_id = relative_parts[-4] # Adjust if group is at a different level

                # Basic validation of folder names
                if not temp_str.endswith('K'):
                    # print(f"Skipping {xvg_file}: Parent directory '{temp_str}' doesn't end with 'K'.")
                    continue
                if not time_str.endswith('ns'):
                    # print(f"Skipping {xvg_file}: Parent directory '{time_str}' doesn't end with 'ns'.")
                    continue

                analysis_type = xvg_file.stem # e.g., 'rmsd', 'gyrate'
                temp = temp_str
                time = time_str

                # print(f"Processing: Group={group_id}, Temp={temp}, Time={time}, Analysis={analysis_type}, File={xvg_file.name}")

                data = parse_xvg(xvg_file)
                if data and data['x_data'].size > 0: # Check if parsing was successful and data exists
                    key = (temp, time)
                    # Store under group_id -> analysis_type -> (temp, time)
                    all_data[group_id][analysis_type][key] = data
                # else:
                #     print(f"Skipping {xvg_file} due to parsing error or no data.")

            else:
                # print(f"Skipping {xvg_file}: Path structure not recognized or too shallow relative to root.")
                pass # Or print a warning

        except Exception as e:
            print(f"Error processing path for {xvg_file}: {e}")
            traceback.print_exc()

    # Convert defaultdict back to regular dict for cleaner output (optional)
    return {k: dict(v) for k, v in all_data.items()}


def calculate_statistics(data_collection, group_id, analysis_type):
    """Calculate statistics and return DataFrame"""
    stats_data = []

    for (temp, time), data in data_collection.items():
        y_data = data['y_data']

        # Skip if y_data is empty or malformed
        if y_data is None or y_data.size == 0:
            # print(f"Warning: Empty y_data for {group_id}/{analysis_type}/{temp}/{time}. Skipping stats.")
            continue

        num_series = y_data.shape[1] if y_data.ndim == 2 else 1
        is_multi_series = num_series > 1

        # Ensure labels list matches the number of series
        labels = data['labels']
        if len(labels) != num_series:
             # print(f"Warning: Label count mismatch for {group_id}/{analysis_type}/{temp}/{time}. Generating defaults.")
             # Generate default labels if mismatch or none provided
             if is_multi_series:
                 labels = [f'Series {i+1}' for i in range(num_series)]
             elif data.get('y_label'):
                 labels = [data['y_label']]
             else:
                 labels = ['Y-Value'] # Fallback default

        half_idx = len(data['x_data']) // 2 # Use x_data length for index

        for i in range(num_series):
            # Select the i-th column, handling 1D/2D cases
            current_y_series = y_data[:, i] if is_multi_series else y_data.flatten()

            # Check for sufficient data points and ignore NaNs
            valid_series = current_y_series[~np.isnan(current_y_series)]
            valid_equil_series = current_y_series[half_idx:][~np.isnan(current_y_series[half_idx:])]

            if valid_series.size == 0:
                # print(f"Warning: No valid data points in series {i} for {group_id}/{analysis_type}/{temp}/{time}. Skipping.")
                continue # Skip this series if all values are NaN or empty

            mean = np.mean(valid_series)
            std = np.std(valid_series)

            if valid_equil_series.size > 0:
                equil_mean = np.mean(valid_equil_series)
                equil_std = np.std(valid_equil_series)
            else: # Handle case where second half is all NaN or too short
                equil_mean = np.nan
                equil_std = np.nan

            label = labels[i]

            stats_data.append({
                'Group': group_id,
                'Analysis': analysis_type,
                'Temperature': temp,
                'Simulation Time': time,
                'Data Series': label, # Renamed from 'Data Type'
                'Mean': mean,
                'Std Dev': std,
                'Equilibrated Mean': equil_mean,
                'Equilibrated Std Dev': equil_std
            })

    return pd.DataFrame(stats_data)


def plot_comparison_by_temperature(data_collection, output_dir, group_id, analysis_type, time_point=None):
    """Compare data across different temperatures for a specific group and analysis."""
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Group data by simulation time first
    time_groups = defaultdict(list)
    for (temp, time), data in data_collection.items():
        time_groups[time].append((temp, data))

    # Sort time groups naturally (e.g., 1ns, 10ns, 50ns)
    sorted_times = sorted(time_groups.keys(), key=lambda t: int(t.replace('ns', '')))

    for time in sorted_times:
        if time_point and time != time_point:
            continue # Skip if a specific time is requested and this isn't it

        temp_data_list = sorted(time_groups[time], key=lambda item: int(item[0].replace('K',''))) # Sort by temp

        if not temp_data_list: continue # Skip if no data for this time point

        plt.figure(figsize=(12, 7))
        # Use data from the first item for labels, assuming consistency
        first_temp, first_data = temp_data_list[0]

        # Check if the first dataset has multiple y-series
        num_series = first_data['y_data'].shape[1] if first_data['y_data'].ndim == 2 else 1
        is_multi_series = num_series > 1

        plot_successful = False
        for temp, data in temp_data_list:
            x_data = data['x_data']
            y_data = data['y_data']

            if x_data.size == 0 or y_data.size == 0: continue # Skip empty datasets

            current_num_series = y_data.shape[1] if y_data.ndim == 2 else 1
            labels = data['labels']
             # Ensure labels list matches the number of series for *this* specific data
            if len(labels) != current_num_series:
                if current_num_series > 1:
                    labels = [f'Series {i+1}' for i in range(current_num_series)]
                elif data.get('y_label'):
                    labels = [data['y_label']]
                else:
                    labels = ['Y-Value']


            if current_num_series > 1: # Multi-series plot
                for i in range(current_num_series):
                    series_label = labels[i]
                    plt.plot(x_data, y_data[:, i], label=f"{temp} - {series_label}", alpha=0.8)
                    plot_successful = True
            elif current_num_series == 1: # Single series plot
                 plt.plot(x_data, y_data.flatten(), label=temp, alpha=0.8)
                 plot_successful = True

        if not plot_successful:
            plt.close() # Close the figure if nothing was plotted
            # print(f"Skipping plot {group_id}_{analysis_type}_temp_comparison_{time}.pdf: No data to plot.")
            continue

        # Construct title and labels
        analysis_title = analysis_type.replace('_', ' ').title()
        plot_title = f"{group_id}: {analysis_title} vs Temperature ({time})"
        plt.title(plot_title)
        plt.xlabel(first_data['x_label'] or "Time (ps)") # Use ps as common GROMACS unit
        plt.ylabel(first_data['y_label'] or analysis_title) # Use analysis type if no y-label
        plt.legend(loc='best', fontsize='small')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()

        output_filename = f"{group_id}_{analysis_type}_temp_comparison_{time}.pdf"
        output_path = output_dir / output_filename
        try:
            plt.savefig(output_path, format='pdf', bbox_inches='tight')
            print(f"Temperature comparison plot saved: {output_path}")
        except Exception as e:
            print(f"Error saving plot {output_path}: {e}")
        plt.close()


def plot_comparison_by_time(data_collection, output_dir, group_id, analysis_type, temperature=None):
    """Compare data across different simulation times for a specific group, analysis, and temperature."""
    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Group data by temperature first
    temp_groups = defaultdict(list)
    for (temp, time), data in data_collection.items():
        temp_groups[temp].append((time, data))

    # Sort temp groups naturally (e.g., 280K, 300K, 320K)
    sorted_temps = sorted(temp_groups.keys(), key=lambda t: int(t.replace('K', '')))

    for temp in sorted_temps:
        if temperature and temp != temperature:
            continue # Skip if a specific temperature is requested and this isn't it

        time_data_list = sorted(temp_groups[temp], key=lambda item: int(item[0].replace('ns',''))) # Sort by time

        if not time_data_list: continue # Skip if no data for this temperature
        plt.figure(figsize=(12, 7))
        # Use data from the first item for labels, assuming consistency
        first_time, first_data = time_data_list[0]

        # Check if the first dataset has multiple y-series
        num_series = first_data['y_data'].shape[1] if first_data['y_data'].ndim == 2 else 1
        is_multi_series = num_series > 1

        plot_successful = False
        for time, data in time_data_list:
            x_data = data['x_data']
            y_data = data['y_data']

            if x_data.size == 0 or y_data.size == 0: continue # Skip empty datasets

            current_num_series = y_data.shape[1] if y_data.ndim == 2 else 1
            labels = data['labels']
             # Ensure labels list matches the number of series for *this* specific data
            if len(labels) != current_num_series:
                if current_num_series > 1:
                    labels = [f'Series {i+1}' for i in range(current_num_series)]
                elif data.get('y_label'):
                    labels = [data['y_label']]
                else:
                    labels = ['Y-Value'] # Fallback default

            if current_num_series > 1: # Multi-series plot
                for i in range(current_num_series):
                    series_label = labels[i]
                    plt.plot(x_data, y_data[:, i], label=f"{time} - {series_label}", alpha=0.8)
                    plot_successful = True
            elif current_num_series == 1: # Single series plot
                 plt.plot(x_data, y_data.flatten(), label=time, alpha=0.8)
                 plot_successful = True

        if not plot_successful:
            plt.close() # Close the figure if nothing was plotted
            # print(f"Skipping plot {group_id}_{analysis_type}_time_comparison_{temp}.pdf: No data to plot.")
            continue

        # Construct title and labels
        analysis_title = analysis_type.replace('_', ' ').title()
        plot_title = f"{group_id}: {analysis_title} vs Time ({temp})"
        plt.title(plot_title)
        plt.xlabel(first_data['x_label'] or "Time (ps)") # Use ps as common GROMACS unit
        plt.ylabel(first_data['y_label'] or analysis_title) # Use analysis type if no y-label
        plt.legend(loc='best', fontsize='small')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()

        output_filename = f"{group_id}_{analysis_type}_time_comparison_{temp}.pdf"
        output_path = output_dir / output_filename
        try:
            plt.savefig(output_path, format='pdf', bbox_inches='tight')
            print(f"Time comparison plot saved: {output_path}")
        except Exception as e:
            print(f"Error saving plot {output_path}: {e}")
        plt.close()


def plot_bar_statistics(stats_df, output_dir, group_id, analysis_type):
    """Create bar plots and heatmaps for statistical data for a specific group and analysis."""
    if stats_df.empty:
        print(f"Skipping statistics plots for {group_id}/{analysis_type}: No data.")
        return

    output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    sns.set_style("whitegrid")

    # --- Plotting Logic ---
    # We might have multiple data series (e.g., from multi-column XVG or different labels)
    # Plot each 'Data Series' separately or consider how to combine/facet them.

    # Option 1: Plot each 'Data Series' on its own bar chart/heatmap
    for series_label, series_df in stats_df.groupby('Data Series'):
        if series_df.empty:
            continue

        # Prepare data for plotting (handle potential missing temp/time combos)
        try:
            pivot_data_mean = series_df.pivot_table(
                index='Temperature', columns='Simulation Time', values='Equilibrated Mean'
            )
            pivot_data_std = series_df.pivot_table(
                index='Temperature', columns='Simulation Time', values='Equilibrated Std Dev'
            )

            # Sort columns (time) and index (temp) naturally
            pivot_data_mean = pivot_data_mean.sort_index(key=lambda x: x.str.replace('K','').astype(int))
            pivot_data_mean = pivot_data_mean.reindex(sorted(pivot_data_mean.columns, key=lambda x: int(x.replace('ns',''))), axis=1)

            pivot_data_std = pivot_data_std.sort_index(key=lambda x: x.str.replace('K','').astype(int))
            pivot_data_std = pivot_data_std.reindex(sorted(pivot_data_std.columns, key=lambda x: int(x.replace('ns',''))), axis=1)

        except Exception as e:
            print(f"Error pivoting data for {group_id}/{analysis_type}/{series_label}: {e}")
            continue # Skip this series if pivoting fails


        # === Bar Plot (Equilibrated Mean with Std Dev Error Bars) ===
        if not pivot_data_mean.empty:
            fig, ax = plt.subplots(figsize=(12, 8))
            try:
                # Plot bars with error bars (using equilibrated std dev)
                pivot_data_mean.plot(kind='bar', yerr=pivot_data_std, ax=ax,
                                     capsize=4, width=0.8,
                                     cmap='viridis') # Use a colormap

                analysis_title = analysis_type.replace('_', ' ').title()
                series_title_part = f" ({series_label})" if stats_df['Data Series'].nunique() > 1 else ""
                plt.title(f"{group_id}: {analysis_title}{series_title_part} - Equilibrated Mean Comparison")

                plt.ylabel(f"Mean Value ({stats_df['Data Series'].iloc[0]})") # Use series label in ylabel
                plt.xlabel("Temperature")
                plt.xticks(rotation=45, ha='right')
                plt.grid(True, linestyle='--', alpha=0.7, axis='y')
                plt.legend(title="Sim Time")
                plt.tight_layout()

                safe_series_label = series_label.replace('/', '_').replace(' ','_') # Sanitize label for filename
                output_filename = f"{group_id}_{analysis_type}_{safe_series_label}_mean_comparison.pdf"
                output_path = output_dir / output_filename
                plt.savefig(output_path, format='pdf', bbox_inches='tight')
                print(f"Statistical bar plot saved: {output_path}")
            except Exception as e:
                print(f"Error creating/saving bar plot for {group_id}/{analysis_type}/{series_label}: {e}")
            finally:
                plt.close(fig) # Ensure figure is closed

        # === Heatmap (Equilibrated Mean) ===
        if not pivot_data_mean.empty and pivot_data_mean.size > 1: # Heatmap needs >1 cell
            plt.figure(figsize=(10, 8))
            try:
                sns.heatmap(pivot_data_mean, annot=True, cmap='viridis', fmt='.3f', linewidths=.5, cbar_kws={'label': 'Equilibrated Mean'})
                analysis_title = analysis_type.replace('_', ' ').title()
                series_title_part = f" ({series_label})" if stats_df['Data Series'].nunique() > 1 else ""
                plt.title(f"{group_id}: {analysis_title}{series_title_part} - Equilibrated Mean Heatmap")
                plt.tight_layout()

                safe_series_label = series_label.replace('/', '_').replace(' ','_') # Sanitize label for filename
                output_filename = f"{group_id}_{analysis_type}_{safe_series_label}_heatmap.pdf"
                output_path = output_dir / output_filename
                plt.savefig(output_path, format='pdf', bbox_inches='tight')
                print(f"Statistical heatmap saved: {output_path}")
            except Exception as e:
                 print(f"Error creating/saving heatmap for {group_id}/{analysis_type}/{series_label}: {e}")
            finally:
                plt.close() # Ensure figure is closed
        # else:
        #     if pivot_data_mean.size <= 1:
        #          print(f"Skipping heatmap for {group_id}/{analysis_type}/{series_label}: Not enough data points for a heatmap.")


def main():
    parser = argparse.ArgumentParser(description='GROMACS XVG Comparison Analysis Tool')
    parser.add_argument('-i','--input', required=True, help='Root directory containing simulation results (e.g., ./chains)')
    parser.add_argument('-o', '--output', required=True, help='Output directory for plots and statistics')
    parser.add_argument('-t', '--temp', help='Specify temperature for time comparison plots (e.g., 300K)')
    parser.add_argument('-s', '--time', help='Specify simulation time for temperature comparison plots (e.g., 50ns)')
    args = parser.parse_args()

    root_input_dir = Path(args.input)
    if not root_input_dir.is_dir():
        print(f"Error: Input directory '{args.input}' not found or is not a directory.")
        sys.exit(1)

    # Collect all data recursively, grouped by the directory structure
    # all_data structure: {group_id: {analysis_type: {(temp, time): data_dict}}}
    all_data = collect_all_data(root_input_dir)

    if not all_data:
        print(f"No valid XVG data found in subdirectories of {root_input_dir}. Exiting.")
        return # Exit if no data collected

    root_output_dir = Path(args.output)
    os.makedirs(root_output_dir, exist_ok=True)
    print(f"\nStarting analysis. Results will be saved in: {root_output_dir}")

    # Iterate through each group found (e.g., 'ChainD')
    for group_id, group_data in all_data.items():
        print(f"\nProcessing Group: {group_id}")
        group_output_dir = root_output_dir / group_id
        os.makedirs(group_output_dir, exist_ok=True)

        # Iterate through each analysis type within the group (e.g., 'rmsd', 'gyrate')
        for analysis_type, data_collection in group_data.items():
            print(f"  Analysis Type: {analysis_type}")

            if not data_collection:
                print(f"    No data found for {analysis_type}. Skipping.")
                continue

            # Create a specific output directory for this analysis type within the group
            analysis_output_dir = group_output_dir / analysis_type
            os.makedirs(analysis_output_dir, exist_ok=True)

            # 1. Calculate Statistics
            stats_df = calculate_statistics(data_collection, group_id, analysis_type)
            if not stats_df.empty:
                csv_filename = f"{group_id}_{analysis_type}_statistics.csv"
                csv_path = analysis_output_dir / csv_filename
                try:
                    stats_df.to_csv(csv_path, index=False)
                    print(f"    Statistics saved to: {csv_path}")
                except Exception as e:
                    print(f"    Error saving statistics CSV {csv_path}: {e}")
            else:
                 print(f"    No statistics generated for {group_id}/{analysis_type}.")


            # 2. Plot Comparison by Temperature
            print(f"    Generating temperature comparison plots...")
            plot_comparison_by_temperature(data_collection, analysis_output_dir, group_id, analysis_type, args.time)

            # 3. Plot Comparison by Time
            print(f"    Generating time comparison plots...")
            plot_comparison_by_time(data_collection, analysis_output_dir, group_id, analysis_type, args.temp)

            # 4. Plot Bar Statistics and Heatmap
            if not stats_df.empty:
                print(f"    Generating statistics plots...")
                plot_bar_statistics(stats_df, analysis_output_dir, group_id, analysis_type)
            else:
                print(f"    Skipping statistics plots for {group_id}/{analysis_type} (no stats data).")


    print("\nAnalysis complete!")

if __name__ == "__main__":
    # Example usage from command line:
    # python your_script_name.py -i /Users/chen_yiru/Desktop/CMML_MD_result/MD_plots_xvg.txt/chains -o ./analysis_output -t 300K -s 50ns
    main()