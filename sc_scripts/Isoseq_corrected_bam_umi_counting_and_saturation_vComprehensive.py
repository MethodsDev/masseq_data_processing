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
    # Create the argument parser.
    parser = argparse.ArgumentParser(
        description="Process a BAM file, count UMIs, and plot sequencing saturation."
    )
    
    # Add arguments for the file paths.
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
    
    bam_file = args.bam_file
    tsv_output = args.tsv_output
    bcstats_file = args.bcstats_file
    
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            try:
                cell_barcode = read.get_tag("CB")
                umi = read.get_tag("XM")  # UMI from "XM" tag.
            except KeyError:
                continue
            read_counts[cell_barcode][umi] += 1
    
    # Write UMI counts to TSV.
    with open(tsv_output, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["cell_barcode", "UMI", "Count"])
        for cell, umi_dict in read_counts.items():
            for umi, count in umi_dict.items():
                writer.writerow([cell, umi, count])
    
    # Compute total reads from all cells (all reads in BAM).
    total_reads_all = sum(sum(umi_dict.values()) for umi_dict in read_counts.values())
    
    # -------------------------------
    # Part 2: Filter Real Cells and Compute Saturation Metrics
    # -------------------------------
    
    bcstats_df = pd.read_csv(bcstats_file, sep="\t")
    # Filter for real cell barcodes (where RealCell column equals "cell")
    real_barcodes = set(bcstats_df.loc[bcstats_df["RealCell"] == "cell", "#BarcodeSequence"])
    
    reads_per_cell = []      # Total reads per filtered (real) cell.
    saturation_indices = []  # Per-cell saturation index.
    unique_umis_per_cell = []  # For global saturation calculation.
    
    for cell, umi_dict in read_counts.items():
        if cell in real_barcodes:
            total_reads = sum(umi_dict.values())
            unique_umis = len(umi_dict)
            unique_umis_per_cell.append(unique_umis)
            saturation = 1 - (unique_umis / total_reads) if total_reads > 0 else 0
            reads_per_cell.append(total_reads)
            saturation_indices.append(saturation)
    
    mean_reads = np.mean(reads_per_cell) if reads_per_cell else 0
    mean_saturation = np.mean(saturation_indices) if saturation_indices else 0
    
    global_total_reads = sum(reads_per_cell)
    global_total_unique = sum(unique_umis_per_cell)
    global_saturation = 1 - (global_total_unique / global_total_reads) if global_total_reads > 0 else 0
    
    filtered_cell_count = len(reads_per_cell)
    
    # Calculate the fraction of reads going to filtered (real) cells.
    fraction_filtered = global_total_reads / total_reads_all if total_reads_all > 0 else 0
    
    print(f"Average reads per real cell: {mean_reads:.1f}")
    print(f"Total filtered (real) cells: {filtered_cell_count}")
    print(f"Fraction of reads going to filtered cells: {fraction_filtered:.3f}")
    
    # -------------------------------
    # Part 3: Model Fitting and Side-by-Side Plotting
    # -------------------------------
    
    x_data = np.array(reads_per_cell)
    y_data = np.array(saturation_indices)
    
    # Fit the Michaelis–Menten model.
    p0 = [max(y_data), np.mean(x_data)]
    best_fit_params, covar = curve_fit(MM, x_data, y_data, p0=p0)
    Vmax, Km = best_fit_params
    sigma = np.sqrt(np.diag(covar))
    
    # Compute R².
    residuals = y_data - MM(x_data, *best_fit_params)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # Generate a high-resolution fitted curve.
    x_fit = np.linspace(np.min(x_data), np.max(x_data), 1000)
    y_fit = MM(x_fit, *best_fit_params)
    
    # --- Refinement: Trim fitted curve up to 90% of Vmax.
    y_max_val = Vmax  # Vmax from the model.
    y_target = y_max_val * 0.9
    x_coef_indices = np.where(y_fit >= y_target)[0]
    if len(x_coef_indices) > 0:
        cutoff_index = x_coef_indices[0]
        x_fit_trim = x_fit[:cutoff_index]
        y_fit_trim = y_fit[:cutoff_index]
    else:
        x_fit_trim = x_fit
        y_fit_trim = y_fit
    
    # --- Additional Calculation: Extra Reads Needed to Reach 10,000 Reads per Cell.
    target_reads = 10000
    extra_per_cell = max(0, target_reads - mean_reads)
    extra_reads_filtered = extra_per_cell * filtered_cell_count
    # Adjust for the fact that filtered cells get only a fraction of total reads.
    total_extra_reads_required = extra_reads_filtered / fraction_filtered if fraction_filtered > 0 else 0
    
    # -------------------------------
    # Create Side-by-Side Subplots (Linear and Log Scale)
    # -------------------------------
    
    fig, (ax_lin, ax_log) = plt.subplots(ncols=2, figsize=(16, 8))
    
    def plot_saturation(ax, x_scale="linear"):
        # Plot the trimmed fitted MM model (with thick black lines).
        ax.plot(x_fit_trim, MM(x_fit_trim, *best_fit_params), color='k', linewidth=3, label=f"Vmax = {Vmax:.3f}")
        ax.plot(x_fit_trim, MM(x_fit_trim, *best_fit_params), color='k', linewidth=3, label=f"Km = {Km:.1f}")
        
        # Confidence intervals.
        bound_upper = MM(x_fit_trim, Vmax + sigma[0], Km + sigma[1])
        bound_lower = MM(x_fit_trim, Vmax - sigma[0], Km - sigma[1])
        ax.fill_between(x_fit_trim, bound_lower, bound_upper, color="gray", alpha=0.4)
        
        # Plot raw data.
        ax.scatter(x_data, y_data, s=10, c="red", label="Real Cells")
        
        # Plot vertical/horizontal lines for mean values.
        ax.axvline(mean_reads, color='red', linestyle='--', label=f"Mean reads per cell: {mean_reads:.1f}")
        ax.axhline(mean_saturation, color='green', linestyle='--', label=f"Avg saturation index: {mean_saturation:.3f}")
        
        # Plot vertical line at the knee point.
        if len(x_coef_indices) > 0:
            knee_x = x_fit[x_coef_indices[0]]
            ax.axvline(knee_x, color='purple', linestyle='--', label=f"Knee point: {knee_x:.0f} reads")
        else:
            knee_x = np.nan
        
        # Annotations.
        ax.text(0.05, 0.95, f"Global Saturation: {global_saturation:.3f}", transform=ax.transAxes,
                fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
        ax.text(0.05, 0.90, f"Mean Reads per Cell: {mean_reads:.1f}", transform=ax.transAxes,
                fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
        ax.text(0.05, 0.85, f"Filtered Cells: {filtered_cell_count}", transform=ax.transAxes,
                fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
        ax.text(0.05, 0.80, f"Vmax = {Vmax:.3f} ± {sigma[0]:.3f}\nKm = {Km:.1f} ± {sigma[1]:.1f}\nR² = {r_squared:.3f}",
                transform=ax.transAxes, fontsize=12, verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
        ax.text(0.05, 0.75, f"Extra reads needed per filtered cell: {extra_per_cell:.0f}\n"
                           f"Total extra for filtered cells: {extra_reads_filtered:.0f}\n"
                           f"Total extra sequencing required: {total_extra_reads_required:.0f}",
                transform=ax.transAxes, fontsize=12, verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
        ax.text(0.05, 0.70, f"Fraction filtered: {fraction_filtered*100:.1f}%", transform=ax.transAxes,
                fontsize=12, verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='gray'))
        
        ax.set_xlabel("Reads per cell")
        ax.set_ylabel("Sequencing Saturation Index")
        ax.set_title(f"X-axis {x_scale.capitalize()} Scale")
        ax.legend()
        ax.grid(True)
        ax.set_xscale(x_scale)
    
    # Plot on left (linear scale).
    plot_saturation(ax_lin, x_scale="linear")
    
    # Plot on right (log scale).
    plot_saturation(ax_log, x_scale="log")
    
    plt.suptitle("Sequencing Saturation for Real Cells", fontsize=16)
    
    # Save and display the figure.
    plt.savefig(args.saturation_index_plot)
    plt.show()

if __name__ == "__main__":
    main()
