#!/usr/bin/env python3
# Roger Volden

'''
Meant to plot the number of reads per cell barcode,
can also plot umis per cell

Usage:
    python3 plot_knees.py \
            -b skera.fltnc.bam
            -o plots/skera
            -t {cbc,umi,both}

    OR

    python3 plot_knees.py \
            -s bcstats.tsv \
            -o plots/skera
'''

import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from palette import Palette as pbp

def parse_args():
    parser = argparse.ArgumentParser()
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        '--tsv', '-s', type=str,
        help='Output tsv file from bcstats (can be gzipped)'
    )
    inputs.add_argument(
        '--bam', '-b', type=str,
        help='Input bam'
    )
    parser.add_argument(
        '--output', '-o', type=str, required=True, help='Output png prefix'
    )
    parser.add_argument(
        '--type', '-t', type=str,
        choices=['cbc', 'umi', 'both'], default='cbc',
        help='Type of plot to output (default both)'
    )
    parser.add_argument(
        '--cumulative', '-c', action='store_true', default=False,
        help='Use to generate cumulative counts [False]'
    )
    parser.add_argument(
        '--max_cells', '-m', default=-1, type=int,
        help='Force an x axis maximum instead of the mpl default'
    )
    return parser.parse_args()

def read_counts(counts):
    '''
    Reads in the dedup csv file and keeps track
    of how many reads per cell barcode.

    read_counts = {'cell_bc_1': [read_count, {'molecule/1', ...}], ...}
    '''
    read_counts = {}
    pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(counts, 'rb', check_sq=False)
    for read in tqdm(bam.fetch(until_eof=True), desc='Reading bam'):
        read_id = read.query_name
        try:
            cbc, umi = read.get_tag('XC'), read.get_tag('XM')
        except KeyError:
            cbc, umi = read.get_tag('CB'), read.get_tag('XM')
        if not cbc or not umi:
            continue
        if cbc not in read_counts:
            read_counts[cbc] = [1, set([umi])]
        else:
            read_counts[cbc][1].add(umi)
            read_counts[cbc][0] += 1
    bam.close()
    return read_counts

def read_tsv(tsv, max_cells):
    zipped, first = False, True
    if tsv.endswith('.gz'):
        import gzip
        fh = gzip.open(tsv, 'rb')
        zipped = True
    else:
        fh = open(tsv, 'r')

    counts, last_real = [], -1
    for line in fh:
        if first:
            first = False
            continue
        if zipped:
            line = line.decode()
        _, nreads, rank, numis, real, _ = line.rstrip().split('\t')
        nreads, rank, numis = int(nreads), int(rank), int(numis)
        if rank >= max_cells and max_cells > 0:
            break
        if real == 'cell':
            last_real = rank
        counts.append(numis)

    fh.close()
    counts.sort(reverse=True)
    return counts, last_real

def plot_from_tsv(counts, ncells, args):
    max_cells, output = args.max_cells, args.output

    plt.figure(figsize=(5, 3))
    # plt.style.use('clean')
    c = plt.axes([0.125, 0.125, 0.8, 0.8])

    colors = [
        pbp.light_orange, pbp.light_pink,  pbp.light_purple, pbp.light_blue,
        pbp.light_teal,   pbp.light_green, pbp.orange,       pbp.pink,
        pbp.purple,       pbp.blue,        pbp.teal,         pbp.green,
        pbp.dark_orange,  pbp.dark_blue,   pbp.dark_green,   pbp.dark_grey
    ]

    # x, y = list(zip(*counts))
    # print(x, y)
    c.plot(range(len(counts)), counts, lw=1, color='grey', zorder=10)
    c.plot(range(ncells), counts[:ncells], lw=1.5, color=pbp.pink, zorder=11)

    c.set_xlabel(r'Cell # (log$_{10}$)')
    c.set_xscale('log')
    c.set_yscale('log')
    c.set_ylabel(r'log$_{10}$(# of UMIs)')
    c.set_title('UMIs per cell')

    output += '.knee.png'
    plt.savefig(output, dpi=600)

def plot_rpc(counts, args):
    max_cells, output = args.max_cells, args.output
    cumulative = args.cumulative

    plt.figure(figsize=(5, 3))
    # plt.style.use('clean')
    c = plt.axes([0.125, 0.125, 0.8, 0.8])

    colors = [
        pbp.light_orange, pbp.light_pink,  pbp.light_purple, pbp.light_blue,
        pbp.light_teal,   pbp.light_green, pbp.orange,       pbp.pink,
        pbp.purple,       pbp.blue,        pbp.teal,         pbp.green,
        pbp.dark_orange,  pbp.dark_blue,   pbp.dark_green,   pbp.dark_grey
    ]

    if cumulative:
        non_c = sorted([x[0] for x in list(counts.values())], reverse=True)
        y = np.cumsum(non_c[:max_cells])
        c.plot(range(len(y)), np.log10(y), lw=1, color='black')
        c.set_xlabel('Cell #')
        c.set_ylabel(r'log$_{10}$(Cumulative # of reads)')
        c.set_title('Cumulative reads per cell')
    else:
        y = sorted([x[0] for x in list(counts.values())], reverse=True)
        umis = sorted([len(x[1]) for x in list(counts.values())], reverse=True)
        cutoff = np.percentile(umis, 99) * 10
        real_cells = [x for x in umis if x >= cutoff]
        up_to_cutoff = len(real_cells)
        print(f'UMI cutoff: {cutoff}')
        print(f'Cells: {up_to_cutoff}')
        print(f'Median UMIs per cell: {np.median(real_cells)}')
        print(f'Total reads up to cutoff: {sum(y[:up_to_cutoff])}')
        print(f'%-age of reads in cells: {sum(y[:up_to_cutoff])/sum(umis)*100:.2f}')
        if max_cells < 0:
            max_cells = len(umis)
        c.plot(
            range(max_cells),
            umis[:max_cells], 
            lw=1, color='grey', zorder=10
        )
        c.plot(
            range(up_to_cutoff),
            umis[:up_to_cutoff], 
            lw=1.5, color=pbp.pink, zorder=11
        )
        c.set_xlabel(r'Cell # (log$_{10}$)')
        c.set_xscale('log')
        c.set_yscale('log')
        c.set_ylabel(r'log$_{10}$(# of UMIs)')
        c.set_title('UMIs per cell')

    output += '.rpc.png'
    plt.savefig(output, dpi=600)

def plot_upc(counts, args):
    max_cells, output = args.max_cells, args.output
    cumulative = args.cumulative

    plt.figure(figsize=(8, 8))
    # plt.style.use('clean')
    c = plt.axes([0.125, 0.125, 0.8, 0.8])

    y = sorted([len(x[1]) for x in list(counts.values())], reverse=True)
    print(f'Total UMIs up to cutoff: {sum(y[:max_cells])}')
    if max_cells < 0:
        c.plot(y, lw=1, color='black')
    else:
        c.plot(
            y[:max_cells], 
            lw=1, color='black'
        )

    c.set_xscale("log")
    c.set_yscale("log")
    c.set_xlabel('Cell #')
    c.set_ylabel('# of UMIs')
    c.set_title('UMIs per cell')

    output += '.upc.png'
    plt.savefig(output)

def main(args):
    if args.bam:
        counts = read_counts(args.bam)
        if args.type == 'cbc':
            plot_rpc(counts, args)
        elif args.type == 'umi':
            plot_upc(counts, args)
        else:
            plot_rpc(counts, args)
            plot_upc(counts, args)
    else:
        counts, ncells = read_tsv(args.tsv, args.max_cells)
        plot_from_tsv(counts, ncells, args)

    # print('#barcode\treads\tumis')
    # for cbc, countlist in tqdm(counts.items()):
    #     print(f'{cbc}\t{countlist[0]}\t{len(countlist[1])}')
    # exit()


if __name__ == '__main__':
    args = parse_args()
    main(args)
