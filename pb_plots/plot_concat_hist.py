#!/usr/bin/env python3
# Roger Volden

'''
Takes a skera output and plots a histogram of the number of molecules
in each concatemer.

Usage:
    python3 plot_concat_hist.py \
        -s skera.read_lengths.csv \
        -o output.png

Old usage:
    python3 plot_concat_hist.py \
        -b skera.out.bam \
        -c len_cutoff [100] \
        -o output.png
'''

import argparse
import pysam
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
import numpy as np
# from palette import Palette as pbp

def parse_args():
    parser = argparse.ArgumentParser()
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        '--csv', '-s', type=str,
        help='Read length csv file from skera (can be gzipped)'
    )
    inputs.add_argument(
        '--bam', '-b', type=str,
        help='Bam file output by Skera'
    )
    parser.add_argument(
        '--arraysize', '-a', type=int, default=15,
        help='Array size [15]'
    )
    parser.add_argument(
        '--cutoff', '-c', type=int, default=100,
        help='Minimum length cutoff [100]'
    )
    parser.add_argument(
        '--output', '-o', type=str, default='skera_concat_hist.png',
        help='Output png file [skera_concat_hist.png]'
    )
    return parser.parse_args()

def read_csv(csv):
    '''
    Reads through the csv file (can be zipped) and only keeps
    the concatenation factor
    count_dict = {zmw: 8, ...}
    '''
    zipped, first = False, True
    if csv.endswith('.gz'):
        import gzip
        fh = gzip.open(csv, 'rb')
        zipped = True
    else:
        fh = open(csv, 'r')

    count_dict = {}
 
    for line in fh:
        if first:
            first = False
            continue
        if zipped:
            line = line.decode()
        zmw, hifi_rl, d_rl, concat = [int(x) for x in line.rstrip().split(',')]

        if zmw in count_dict:
            continue
        count_dict[zmw] = concat

    fh.close()
    return count_dict

def make_count_dict(bam, cutoff):
    '''
    Reads through the Skera output to make a dictionary of header 
    to how many reads it got split into

    count_dict = {whatever/ccs: 8, ...}
    '''
    count_dict = {}
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'rb', check_sq=False, require_index=False) as bam_file:
        for read in tqdm(bam_file.fetch(until_eof=True), desc='Reading skera bam'):
            read = read.tostring(bam_file).split('\t')
            name, s_len = read[0], len(read[9])
            if s_len < cutoff: continue
            if name not in count_dict:
                count_dict[name] = 0
            count_dict[name] += 1
    return count_dict

def plot_concat(bins, percents, arraysize, output):
    plt.figure(figsize=(6, 4))
    # plt.style.use('clean')
    h = plt.axes([0.15, 0.125, 0.8, 0.8])
    
    for i in range(len(bins)):
        left = i + 0.5
        box = Rect((left, 0), 1, bins[i], lw=0, fc="blue")
        h.add_patch(box)
        h.text(i+1, bins[i], str(percents[i]) + '%',
            ha='center', va='bottom', fontsize='small',
        )

    h.set_xlim(0.25, arraysize + 0.75)
    h.set_ylim(0, max(bins)*1.1)
    h.set_xlabel('Concatemers per molecule')
    h.set_ylabel('Count')
    h.set_xticks(range(1, arraysize + 1))
    h.set_xticklabels(range(1, arraysize + 1))

    # if not output.endswith('.png'):
    #     output += '.png'
    plt.savefig(output, dpi=600)

def main(args):
    if args.csv:
        count_dict = read_csv(args.csv)
    else:
        count_dict = make_count_dict(args.bam, args.cutoff)
    heights = [0] * args.arraysize
    for molecules in count_dict.values():
        heights[molecules-1] += 1
    total = len(count_dict)
    percents = [round(x/total*100, 2) for x in heights]
    plot_concat(heights, percents, args.arraysize, args.output)

if __name__ == '__main__':
    args = parse_args()
    main(args)
