#!/usr/bin/env python
import argparse
import pysam
from collections import defaultdict
import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Process a BAM file, count UMIs, and plot sequencing saturation.")
    
    # Add arguments for the file paths (non-positional, with flags)
    parser.add_argument('--bam_file', type=str, required=True, help="Input BAM file path")
    parser.add_argument('--tsv_output', type=str, required=True, help="Output TSV file path for UMI counts")
    parser.add_argument('--bcstats_file', type=str, required=True, help="Barcode stats file path with cell metadata")
    parser.add_argument('--saturation_index_plot', type=str, required=True, help="Output PNG file path for the saturation index plot")
    
    return parser.parse_args()

def MM(x, Vmax, Km):
    """Michaelis–Menten model function."""
    return (Vmax * x) / (Km + x)

def main():
    args = parse_arguments()
    
    # -------------------------------
    # Part 1: Process BAM File and Write UMI Counts to TSV
    # -------------------------------
    
    # Dictionary for per-cell UMI counts:
    # Key: cell barcode (from tag "CB")
    # Value: dictionary with key = UMI (from tag "XM") and value = count
    read_counts = defaultdict(lambda: defaultdict(int))
    
    bam_file = args.bam_file          # Input BAM file
    tsv_output = args.tsv_output        # Output TSV file for UMI counts
    bcstats_file = args.bcstats_file    # Barcode stats file with cell metadata
    
    # Process the BAM file (use check_sq=False to ignore missing SQ header lines)
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            try:
                cell_barcode = read.get_tag("CB")
                # Use the "XM" tag for the UMI (instead of "UB")
                umi = read.get_tag("XM")
            except KeyError:
                continue
            read_counts[cell_barcode][umi] += 1
    
    # Write the UMI counts to a TSV file (columns: cell_barcode, UMI, Count)
    with open(tsv_output, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["cell_barcode", "UMI", "Count"])
        for cell, umi_dict in read_counts.items():
            for umi, count in umi_dict.items():
                writer.writerow([cell, umi, count])
    
    # -------------------------------
    # Part 2: Filter Real Cells from bcstats_report.tsv and Compute Saturation
    # -------------------------------
    
    bcstats_df = pd.read_csv(bcstats_file, sep="\t")
    # Filter for "real" cell barcodes (where RealCell column equals "cell")
    real_barcodes = set(bcstats_df.loc[bcstats_df["RealCell"] == "cell", "#BarcodeSequence"])
    
    reads_per_cell = []      # Total reads per cell
    saturation_indices = []  # Per-cell saturation index
    unique_umis_per_cell = []  # For global saturation calculation
    
    for cell, umi_dict in read_counts.items():
        if cell in real_barcodes:
            total_reads = sum(umi_dict.values())
            unique_umis = len(umi_dict)
            unique_umis_per_cell.append(unique_umis)
            saturation = 1 - (unique_umis / total_reads) if total_reads > 0 else 0
            reads_per_cell.append(total_reads)
            saturation_indices.append(saturation)
    
    mean_reads = sum(reads_per_cell) / len(reads_per_cell) if reads_per_cell else 0
    mean_saturation = sum(saturation_indices) / len(saturation_indices) if saturation_indices else 0
    global_total_reads = sum(reads_per_cell)
    global_total_unique = sum(unique_umis_per_cell)
    global_saturation = 1 - (global_total_unique / global_total_reads) if global_total_reads > 0 else 0
    filtered_cell_count = len(reads_per_cell)
    
    print(f"Average reads per real cell: {mean_reads:.1f}")
    print(f"Total filtered (real) cells: {filtered_cell_count}")
    
    # -------------------------------
    # Part 3: Model Fitting and Enhanced Plotting for Real Cells
    # -------------------------------
    
    # Convert lists to numpy arrays for model fitting.
    x_data = np.array(reads_per_cell)
    y_data = np.array(saturation_indices)
    
    # Determine if x_data values are very high and need a log scale for plotting.
    log_threshold = 1e6  # Threshold can be adjusted as needed.
    use_log = np.max(x_data) > log_threshold
    
    # Fit the Michaelis–Menten model using curve_fit.
    p0 = [max(y_data), np.mean(x_data)]
    best_fit_params, covar = curve_fit(MM, x_data, y_data, p0=p0)
    Vmax, Km = best_fit_params
    sigma = np.sqrt(np.diag(covar))
    
    residuals = y_data - MM(x_data, *best_fit_params)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # Generate a high-resolution fitted curve.
    x_fit = np.linspace(np.min(x_data), np.max(x_data), 1000)
    y_fit = MM(x_fit, *best_fit_params)
    
    # --- Additional Refinements ---
    # Trim the fitted curve to only plot up to 90% of Vmax.
    y_max_val = Vmax  # Maximum model value.
    y_target = y_max_val * 0.9
    x_coef_indices = np.where(y_fit >= y_target)[0]
    if len(x_coef_indices) > 0:
        cutoff_index = x_coef_indices[0]
        x_fit_trim = x_fit[:cutoff_index]
        y_fit_trim = y_fit[:cutoff_index]
        # bug:  x_max_val = np.max(x_fit_trim)
        if x_fit_trim.size > 0:
            x_max_val = np.max(x_fit_trim)
        else:
            print("Warning: x_fit_trim is empty, using default value")
            x_max_val = 0  # Choose appropriate default for your use case
    else:
        x_fit_trim = x_fit
        y_fit_trim = y_fit
        if x_fit > 0 :
            x_max_val = np.max(x_fit)
        else:
            x_max_val = 0 
    
    # Compute confidence intervals using parameter uncertainties.
    bound_upper = MM(x_fit_trim, Vmax + sigma[0], Km + sigma[1])
    bound_lower = MM(x_fit_trim, Vmax - sigma[0], Km - sigma[1])
    
    # -------------------------------
    # Plotting Section
    # -------------------------------
    plt.figure(figsize=(8, 6))
    
    # Plot the refined fitted curve as thick black lines with parameter labels.
    plt.plot(x_fit_trim, MM(x_fit_trim, *best_fit_params), color='k', linewidth=3, label=f"Vmax = {Vmax:.3f}")
    plt.plot(x_fit_trim, MM(x_fit_trim, *best_fit_params), color='k', linewidth=3, label=f"Km = {Km:.1f}")
    
    # Plot confidence intervals.
    plt.fill_between(x_fit_trim, bound_lower, bound_upper, color="gray", alpha=0.4)
    
    # Plot raw data (in red).
    plt.scatter(x_data, y_data, s=10, c="red", label="Real Cells")
    
    # Plot vertical and horizontal lines for overall means.
    plt.axvline(mean_reads, color='red', linestyle='--', label=f"Mean reads per cell: {mean_reads:.1f}")
    plt.axhline(mean_saturation, color='green', linestyle='--', label=f"Avg saturation index: {mean_saturation:.3f}")
    
    # Plot a vertical line at the knee point.
    if len(x_coef_indices) > 0:
        knee_x = x_fit[cutoff_index]
        plt.axvline(knee_x, color='purple', linestyle='--', label=f"Knee point: {knee_x:.0f} reads")
    else:
        knee_x = np.nan
    
    # Annotate the plot with global metrics and model parameters.
    plt.text(0.05, 0.95, f"Global Saturation: {global_saturation:.3f}", transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
    plt.text(0.05, 0.90, f"Mean Reads per Cell: {mean_reads:.1f}", transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
    plt.text(0.05, 0.85, f"Filtered Cells: {filtered_cell_count}", transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
    plt.text(0.05, 0.80, f"Vmax = {Vmax:.3f} ± {sigma[0]:.3f}\nKm = {Km:.1f} ± {sigma[1]:.1f}\nR² = {r_squared:.3f}",
             transform=plt.gca().transAxes, fontsize=12, verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
    
    plt.xlabel("Reads per cell")
    plt.ylabel("Sequencing Saturation Index")
    plt.title("Sequencing Saturation for Real Cells")
    
    # If values are very high, use a log scale on the x-axis.
    if use_log:
        plt.xscale("log")
    
    plt.legend()
    plt.grid(True)
    
    # Save and display the plot.
    plt.savefig(args.saturation_index_plot)
    plt.show()

if __name__ == "__main__":
    main()

