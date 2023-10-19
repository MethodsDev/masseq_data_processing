#!/usr/bin/env python3
# Roger Volden

'''
Takes a csv with ligation info and plots it as a heatmap

Usage:
    python3 plot_ligation_heatmap.py \
        -c ligations.csv \
        -o output.png
'''

import argparse
import pysam
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rect
import numpy as np
from palette import Palette as pbp

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--csv', '-c', required=True, type=str,
        help='Csv containing ligation information'
    )
    parser.add_argument(
        '--output', '-o', type=str, default='skera_ligations.png',
        help='Output png file [skera_ligations.png]'
    )
    parser.add_argument(
        '--arraysize', '-a', type=int, default=15,
        help='Array size [15]'
    )
    parser.add_argument(
        '--exclude_nums', '-e', action='store_true', default=False,
        help='Use to remove the numbers from the heatmap'
    )
    return parser.parse_args()

def read_csv(csv, arraysize):
    '''
    Placeholder for csv parsing (update when I have the file)
    I guess I can just... return a 2d array with the counts?
    
    Current tmp format:
    Adapter_1,Adapter_2,Ligations
    0,0,0
    0,1,472
    ...
    '''
    max_idx, counts = 0, np.zeros((arraysize+1, arraysize+1), dtype=int)
    first = True
    with open(csv) as f:
        for line in f:
            line = line.rstrip()
            if not line or first:
                first = False
                continue
            i, j, n = [int(x) for x in line.split(',')]
            counts[j][i] = n
    return counts

def plot_ligations(counts, exclude_nums, output):
    plt.figure(figsize=(6, 5))
    # plt.style.use('clean')
    h = plt.axes([0.1, 0.1, 4.25/6, 4.25/5])
    scale = plt.axes([0.125+(4.25/6), 0.1, 0.06, 4.25/5])

    edge_len = len(counts)

    # normalize counts from 0 to 100
    counts_norm = np.array(counts)
    # print(counts_norm[::-1])
    c_min, c_max = np.min(counts_norm), np.max(counts_norm)
    counts_norm = ((counts_norm / np.max(counts_norm)) * 100).astype(int)
    # print(counts_norm[::-1])

    # color gradient
    # gradient = [pbp.white, pbp.green]
    # gradient = [pbp.white, pbp.light_orange, pbp.blue]
    # gradient = [pbp.white, pbp.green, pbp.orange, pbp.blue]
    # gradient = [pbp.white, pbp.light_purple, pbp.teal, pbp.light_orange, pbp.blue]
    # gradient = ["white", "light_purple", "teal", "light_orange", "blue"]
    # gradient = plt.cm.
    gradient = [
        pbp.white, pbp.light_pink, pbp.light_orange,
        pbp.light_green, pbp.light_teal, pbp.light_blue, pbp.purple
    ]
    R, G, B = [], [], []
    steps = 100 // (len(gradient)-1)
    difference = 101 - (steps*(len(gradient)-1))
    for i in range(len(gradient)-1):
        if i == len(gradient) - 2:
            steps += difference
        for c in np.linspace(gradient[i][0], gradient[i+1][0], steps):
            R.append(c)
        for c in np.linspace(gradient[i][1], gradient[i+1][1], steps):
            G.append(c)
        for c in np.linspace(gradient[i][2], gradient[i+1][2], steps):
            B.append(c)

    # pick the fontsize
    if len(str(c_max)) <= 5:
        fsize = 5
    elif len(str(c_max)) == 6:
        fsize = 4
    else:
        fsize = 3

    # plot the heatmap
    # plots in reversed row order (first row on the bottom)
    for i in range(edge_len): # row
        for j in range(edge_len): # column
            color = (
                R[counts_norm[j][i]],
                G[counts_norm[j][i]],
                B[counts_norm[j][i]]
            )
            box = Rect((i, j), 1, 1, lw=0, zorder=9, fc=color)
            h.add_patch(box)
            if exclude_nums:
                continue
            h.text(
                i+0.5, j+0.5, str(counts[j][i]),
                ha='center', va='center', fontsize=fsize, zorder=10
            )

    # plot out the scale colors
    for i in range(len(R)):
        box = Rect((0, i), 1, 1, lw=0, fc=(R[i], G[i], B[i]))
        scale.add_patch(box)

    # heatmap params
    letters = [chr(x+65) for x in range(edge_len)]
    h.set_xlim(0, edge_len)
    h.set_ylim(0, edge_len)
    h.set_xticks(np.arange(0.5, edge_len+0.5))
    h.set_yticks(np.arange(0.5, edge_len+0.5))
    h.set_xticklabels(letters)
    h.set_yticklabels(letters)
    h.set_xlabel('Adapter 1')
    h.set_ylabel('Adapter 2')
    h.tick_params(
        axis='both', which='both',
        bottom=True, labelbottom=True,
        left=True, labelleft=True,
        right=False, labelright=False,
        top=False, labeltop=False
    )

    # scale params
    scale.set_xlim(0, 1)
    scale.set_ylim(0, 100)
    scale.set_yticks(range(0, 101, 20))
    scale.set_yticklabels([int(x) for x in np.linspace(0, c_max, 6)], fontsize=6)
    scale.tick_params(
        axis='both', which='both',
        bottom=False, labelbottom=False,
        left=False, labelleft=False,
        right=True, labelright=True,
        top=False, labeltop=False
    )

    if not output.endswith('.png'):
        output += '.png'
    plt.savefig(output, dpi=600)

def main(args):
    counts = read_csv(args.csv, args.arraysize)
    plot_ligations(counts, args.exclude_nums, args.output)

if __name__ == '__main__':
    args = parse_args()
    main(args)
