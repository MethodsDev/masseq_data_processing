#!/usr/bin/env python3
# Roger Volden

'''
Takes a skera output and plots a histogram of the original readlen.
The bar in the histogram will be colored according to the distribution
of the concatenation factors.

Need: 
    - concat factor
    - len of original read

Usage:
    python3 plot_readlen_hist.py \
        -s skera.read_lengths.csv \
        -o output.png

'''

import argparse
import sys
import os
import pysam
import pickle
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
import numpy as np
from palette import Palette as pbp

def parse_args():
    parser = argparse.ArgumentParser()
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        '--csv', '-s', type=str,
        help='Read length csv file from skera (can be gzipped)'
    )
    inputs.add_argument(
        '--ccs', '-C', type=str,
        help='Bam file output by CCS'
    )
    parser.add_argument(
        '--arraysize', '-a', type=int, default=15,
        help='Array size [15]'
    )
    parser.add_argument(
        '--xmax', '-x', type=int, default=25000,
        help='Manually set x-axis limit [25000]'
    )
    parser.add_argument(
        '--bam', '-b', type=str,
        help='Bam file output by Skera'
    )
    parser.add_argument(
        '--pickle', '-p', type=str, default='len_to_concat.pickle',
        help='Pickle filename [len_to_concat.pickle]'
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

def check_args(args):
    if args.csv and args.bam:
        print(
            'Need to use either the CSV or (CCS BAM + deconcat BAM)',
            file=sys.stderr
        )
        exit(1)
    if not args.csv and not args.bam:
        print(
            'To not use the CSV, you need to supply a CCS BAM and a deconcat BAM',
            file=sys.stderr
        )
        exit(1)

def read_csv(csv):
    '''
    Reads through the csv file (can be zipped) and keeps
    the hifi read length and the concatenation factor
    '''
    zipped, first = False, True
    if csv.endswith('.gz'):
        import gzip
        fh = gzip.open(csv, 'rb')
        zipped = True
    else:
        fh = open(csv, 'r')

    used_zmws, len_to_concat = set(), []
    
    for line in fh:
        if first:
            first = False
            continue
        if zipped:
            line = line.decode()
        zmw, hifi_rl, d_rl, concat = [int(x) for x in line.rstrip().split(',')]

        if zmw in used_zmws:
            continue
        len_to_concat.append((hifi_rl, concat))
        used_zmws.add(zmw)

    fh.close()
    return len_to_concat

def make_length_dict(bam):
    '''
    Reads through the ccs output to make a dictionary of header 
    to how long the read is

    length_dict = {whatever/ccs: 10000, ...}
    '''
    length_dict = {}
    pysam.set_verbosity(0)
    with pysam.AlignmentFile(bam, 'rb', check_sq=False, require_index=False) as bam_file:
        for read in tqdm(bam_file.fetch(until_eof=True), desc='Reading skera bam'):
            read = read.tostring(bam_file).split('\t')
            name, s_len = read[0], len(read[9])
            length_dict[name] = s_len
    return length_dict

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

def plot_concat(len_to_concat, output, arraysize, xmax):
    plt.figure(figsize=(6, 4))

    h = plt.axes([0.15, 0.125, 0.65, 0.8])
    legend = plt.axes([0.8125, 0.125, 0.1, 0.8], frameon=False)
 
    binsize = 250
    nbins = xmax // binsize
    bins = np.arange(0, xmax, binsize)
    lens = np.array([x[0] for x in len_to_concat])
    clipped = np.clip(lens, 0, xmax-1)
    heights = [0] * nbins

    per_bin_concat = {} # [bin] = [15] <- dist of concat
    for i, idx in enumerate(clipped//binsize):
        heights[idx] += 1
        if idx not in per_bin_concat:
            per_bin_concat[idx] = np.zeros(arraysize, dtype=int)
        per_bin_concat[idx][len_to_concat[i][1]-1] += 1

    colors = [
        pbp.light_orange, pbp.light_pink,  pbp.light_purple, pbp.light_blue,
        pbp.light_teal,   pbp.light_green, pbp.orange,       pbp.pink,
        pbp.purple,       pbp.blue,        pbp.teal,         pbp.green,
        pbp.dark_orange,  pbp.dark_blue,   pbp.dark_green,   pbp.dark_grey
    ][:arraysize]

    for i, height in enumerate(heights):
        left = i * binsize
        try:
            last_top = 0
            for j, concat in enumerate(per_bin_concat[i]):
                box = Rect(
                    (left, last_top), binsize, concat, 
                    lw=0.1, ec='black', fc=colors[j]
                )
                h.add_patch(box)
                last_top += concat
        except KeyError:
            continue

    h.set_xlim(0, xmax + binsize)
    h.set_ylim(0, max(heights)*1.1)
    h.set_xlabel('Read length, bp')
    h.set_ylabel('Number of Reads')
    h.set_xticks(range(0, xmax+binsize, binsize*10))

    # legend
    for idx, c in enumerate(colors):
        bottom = idx+0.2 + (0.1*idx)
        cbox = Rect((0.2, bottom), 1, 1, lw=0, color=c)
        legend.add_patch(cbox)
        legend.text(1.4, bottom+0.5, str(idx+1) + 'x', ha='left', va='center')
    legend.set_xlim(0, 2)
    legend.set_ylim(0, bottom + 1.2)
    legend.tick_params(
        axis='both', which='both',
        bottom=False, labelbottom=False,
        left=False, labelleft=False,
        right=False, labelright=False,
        top=False, labeltop=False
    )

    plt.savefig(output, dpi=600)


def main(args):
    if args.csv:
        len_to_concat = read_csv(args.csv)
    else:
        if os.path.exists(args.pickle):
            len_to_concat = pickle.load(open(args.pickle, 'rb'))
        else:
            len_dict = make_length_dict(args.ccs)
            count_dict = make_count_dict(args.bam, args.cutoff)
            len_to_concat = []
            for name, concat in count_dict.items():
                tup = (len_dict[name], concat)
                len_to_concat.append(tup)
            pickle.dump(len_to_concat, open(args.pickle, 'wb'))
    plot_concat(len_to_concat, args.output, args.arraysize, args.xmax)

if __name__ == '__main__':
    args = parse_args()
    check_args(args)
    main(args)
