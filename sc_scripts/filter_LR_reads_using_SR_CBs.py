#!/usr/bin/env python3

import pysam
import argparse
from datetime import datetime
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def read_barcodes(barcode_file):
    with open(barcode_file) as f:
        return set(line.strip() for line in f if line.strip())

def plot_top_barcodes(barcode_counts, output_image, top_n=20):
    if not barcode_counts:
        print("No barcode counts to plot.")
        return

    top_bcs = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    labels, counts = zip(*top_bcs)

    plt.figure(figsize=(10, 6))
    plt.bar(range(len(labels)), counts, color='steelblue')
    plt.xticks(range(len(labels)), labels, rotation=90)
    plt.xlabel("Barcode")
    plt.ylabel("Read Count")
    plt.title(f"Top {top_n} Barcodes by Read Count")
    plt.tight_layout()
    plt.savefig(output_image)
    plt.close()

def write_barcode_counts_tsv(barcode_counts, output_prefix):
    with open(output_prefix, 'w') as f:
        f.write("Barcode\tReadCount\n")
        for barcode, count in sorted(barcode_counts.items(), key=lambda x: -x[1]):
            f.write(f"{barcode}\t{count}\n")
    print(f"Barcode count table saved to: {output_prefix}")

def filter_bam_by_barcode(input_bam, output_bam, barcode_file, report_file, sample_id=None):
    barcodes = read_barcodes(barcode_file)
    print(f"Loaded {len(barcodes)} barcodes from list.")

    in_bam = pysam.AlignmentFile(input_bam, "rb", check_sq=False)
    out_bam = pysam.AlignmentFile(output_bam, "wb", template=in_bam)

    total_reads = 0
    retained_reads = 0
    found_barcodes = set()
    barcode_counts = defaultdict(int)

    for read in in_bam:
        total_reads += 1
        try:
            cb = read.get_tag("CB")
        except KeyError:
            continue

        if cb in barcodes:
            matched_cb = cb
        else:
            rc_cb = reverse_complement(cb)
            if rc_cb in barcodes:
                matched_cb = rc_cb
            else:
                continue

        out_bam.write(read)
        retained_reads += 1
        found_barcodes.add(matched_cb)
        barcode_counts[matched_cb] += 1

    in_bam.close()
    out_bam.close()

    fraction_retained = retained_reads / total_reads if total_reads else 0

    # Stats on barcode counts
    if barcode_counts:
        max_cb = max(barcode_counts.items(), key=lambda x: x[1])
        min_cb = min(barcode_counts.items(), key=lambda x: x[1])
        avg_reads_per_cb = sum(barcode_counts.values()) / len(barcode_counts)
        median_reads_per_cb = np.median(list(barcode_counts.values()))
    else:
        max_cb = ("NA", 0)
        min_cb = ("NA", 0)
        avg_reads_per_cb = 0.0

    summary_lines = []
    if sample_id:
        summary_lines.append(f"Report for Sample: {sample_id}")
    else:
        summary_lines.append("Report")

    summary_lines.append(f"Generated: {datetime.now().isoformat()}")
    summary_lines += [
        f"Input BAM: {input_bam}",
        f"Output BAM: {output_bam}",
        f"Barcode file: {barcode_file}",
        f"Total reads: {total_reads}",
        f"Reads retained: {retained_reads}",
        f"Fraction retained: {fraction_retained:.4f}",
        f"Barcodes in list: {len(barcodes)}",
        f"Barcodes found in BAM: {len(found_barcodes)}",
        f"Average reads per retained barcode: {avg_reads_per_cb:.2f}",
        f"Median reads per retained barcode: {median_reads_per_cb:.2f}",
        f"Barcode with most reads: {max_cb[0]} ({max_cb[1]} reads)",
        f"Barcode with fewest reads: {min_cb[0]} ({min_cb[1]} reads)",
    ]

    plot_file = output_bam + ".barcode_plot.png"
    plot_top_barcodes(barcode_counts, plot_file)
    summary_lines.append(f"Plot saved as: {plot_file}")

    barcodes_tsv=output_bam + ".barcode_counts.tsv"
    write_barcode_counts_tsv(barcode_counts, barcodes_tsv)

    with open(report_file, 'w') as rf:
        rf.write('\n'.join(summary_lines) + '\n')

    print("\n".join(summary_lines))


def main():
    parser = argparse.ArgumentParser(description="Filter BAM by barcode list (CB tag).")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output_bam", required=True, help="Output filtered BAM file")
    parser.add_argument("-b", "--barcode_file", required=True, help="File with list of barcodes")
    parser.add_argument("-r", "--report_file", required=True, help="Output report text file")
    parser.add_argument("-s", "--sample_id", required=False, help="Sample ID to include in report header")
    args = parser.parse_args()

    filter_bam_by_barcode(
        input_bam=args.input_bam,
        output_bam=args.output_bam,
        barcode_file=args.barcode_file,
        report_file=args.report_file,
        sample_id=args.sample_id
    )

if __name__ == "__main__":
    main()
