#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:52:03 2024

@author: akhorgad
"""

#-----------------------changelog-----------------------#
##1. 03-18-2024: Added plot counts function 
##2. 03-23-2024: Remove logic to read gzipped lima.counts as it's not the default out
##3. 04-24-2024: Renamed idmap colnames
##             - 10x_barcode_id : IsoSeq_primer
##             - pb_barcode : Kinnex_Adapter
##             - sampleid : Sample_ID
## Queued: add help - note on structure of idmap.csv
## Queued: make merging on colnames case insensitive
## Queued: add pysam func to merge bams
##4. 05-06-2024: Added pysam functions to merge and sort
## Queued: add function to bgzip all bams 
##5. 06-28-2024: delete locally copied refine bams for space optimization
##6. 07-27-2024: Removed sorting - was not required
##               Renamed merged bam to merged.unaligned.bam
##               Added func(sanity_idmap) to perform sanity checks on input idmap file


import os
import re
import glob
import time
import pysam
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
pd.options.display.max_columns = None

# Script objective:
###    1.to take read counts from lima demux results to plot sample wise distribution 
###    2.to merge refined cleaned demuxed reads from technical replicates for each sample in a batch
###    3.rename merged reads referrencing barcodes as provided in 
###       Format of idmap.csv : Required columns ['Sample_ID','Kinnex_Adapter','IsoSeq_primer'] - case insensitive, whitespace tolerant 


####-- changelog.6 --####
def sanity_idmap(idmap_csv):
    idmap=pd.read_csv(idmap_csv, sep=',')
    print("Reading idmap as provided...")
    print(idmap)
    idmap.columns = idmap.columns.str.lower()
    idmap.rename(columns=lambda x: x.strip(), inplace=True)
    idmap.rename(columns=lambda x: x.replace(' ', '_'), inplace=True)
    if 'sample_id' not in idmap.columns:
        raise Exception("Sample Id column not found, or incorrect format. \n \
                        Accepted formats: [Sample_ID][sample_id][Sample ID][sample id] - case insensitive, whitespace tolerant")
    if 'kinnex_adapter' not in idmap.columns:
        raise Exception("Kinnex Adapter column not found, or incorrect format. \n \
                        Accepted formats: [Kinnex_Adapter][kinnex_adapter][Kinnex Adapter][kinnex adapter] - case insensitive, whitespace tolerant")
    if 'isoseq_primer' not in idmap.columns:
        raise Exception("Isoseq primer column not found, or incorrect format. \n \
                        Accepted formats: [IsoSeq_primer][isoseq_primer][Isoseq Primer][isoseq primer] - case insensitive, whitespace tolerant")
    idmap = idmap[['sample_id', 'kinnex_adapter', 'isoseq_primer']]

    if idmap.isnull().values.any():
        raise Exception("Cannot accept empty values in idmap.csv columns 'Sample_ID','Kinnex_Adapter','IsoSeq_primer'\n \
                        Empty values found at row index: %s, col index: %s" % np.where(pd.isnull(idmap)))
    idmap['sample_id'] = idmap['sample_id'].apply(lambda x: x.strip())                    
    idmap['sample_id'] = idmap.sample_id.str.replace(' ', '_').str.replace('.', '_')
    print(idmap)
    return (idmap)


def lima_cnt_combine(dirpath, opath, idmap_df):
    # read the sample id to barcode correspondence
    print("Reading idmap post clean-up/sanity...")
    print(idmap_df.head(3))
    # List of your file paths
 #   file_paths = sorted(glob.glob(dirpath)) 
    lima_path=os.path.join(dirpath,"*.lima.counts")
    print("Lima outs path provided...")
    print(lima_path)
    file_paths = sorted(glob.glob(lima_path)) 
    print(file_paths)
    cdf=pd.DataFrame()
    for f in file_paths:
        ids_list = extract_movie_name(f)
        ####-- changelog.2 --####   
        df = pd.read_csv(f, sep='\t')
        # adding movie name and barcodes extracted from filename
        df['movie_name']=ids_list[1]
        df['pbbc']=ids_list[2] 
        df['bc10x']=df.IdxFirstNamed.apply(lambda x: re.match(("(bc\d+)"),x).groups()[0] )
        cdf=pd.concat([cdf, df[['movie_name','pbbc','bc10x','Counts']]])  
    print(cdf.head(5))

    # creating a combined counts file for each dataframe
    cdfm=cdf.merge(idmap_df, how='left', left_on=['bc10x','pbbc'], right_on=['isoseq_primer','kinnex_adapter'])
    cdfm=cdfm.drop(columns=['pbbc', 'bc10x'])
    key = lambda x: (x !='movie_name', x != 'sample_id', x != 'isoseq_primer', x != 'kinnex_adapter')
    cdfm=cdfm[sorted(cdfm, key=key)]
    print("Write lima_counts_by_moviename.tsv to..")
    print(os.path.join(opath,"lima_counts_by_moviename.tsv"))
    cdfm.to_csv(os.path.join(opath,"lima_counts_by_moviename.tsv"),sep='\t', index=False)  

    # aggregate counts for replicates assuming these are technical replicates - prepped using same barcodes across multiple flowcells
    agg_bybcs_df=cdfm.groupby(['sample_id','kinnex_adapter','isoseq_primer'])['counts'].agg('sum')
    agg_bybcs_dff=agg_bybcs_df.reset_index()
    print(agg_bybcs_dff)
    print(' Write out aggregated_lima_counts_by_sample.tsv to...')
    print (os.path.join(opath,"aggregated_lima_counts_by_sample.tsv"))
    agg_bybcs_dff.to_csv(os.path.join(opath,"aggregated_lima_counts_by_sample.tsv"),sep='\t', index=False)  
    return(os.path.join(opath,"aggregated_lima_counts_by_sample.tsv"))


# take file lima.counts path and return list of moviename and kinnex bc
def extract_movie_name(filepath):
    hifi_id=re.search(r'([\w.-]+).lima.lima.counts', os.path.basename(filepath)).group(1)
    movie_name=hifi_id.split('.')[0]
    pb_id=hifi_id.split('.')[1]    
    return ([hifi_id,movie_name,pb_id])


def plot_counts(cntpath, opath, plt_title_str="Read counts distribution by sample"):
    tc=pd.read_csv(cntpath, sep='\t')
    tc.head()
    fig, ax = plt.subplots(figsize=(9,4))
    sns.scatterplot(x="sample_id", y="counts", data=tc.sort_values(by=['kinnex_adapter']),  hue="kinnex_adapter", palette='colorblind').tick_params(axis='x', rotation=90)
    plt.grid(color='grey', linestyle='-', linewidth=0.2)
    ax.set(xlabel="Sample name", ylabel=" Read Counts", title=plt_title_str)
    plt.savefig(os.path.join(opath,"readcounts_by_sample.png"), dpi=300, bbox_inches = "tight")
    pass


def merge_bams(grouped_bams_list, op_bam_id, outpath):
    print(grouped_bams_list)
    print(op_bam_id)
    if os.path.exists(outpath):
        opath=os.path.join(outpath,str(op_bam_id+".merged.unaligned.bam"))
    print(opath)
    pysam.merge("-f","-o", opath, *grouped_bams_list)
    ####-- changelog.5 --#### 
    time.sleep(60) 
    for bam in grouped_bams_list: 
        if os.path.isfile(bam):
            os.remove(bam)
        else:
            print("Could not delete {} , file not found. ".format(bam))
    return   


def group_bams(lima_counts_by_moviename_path, bampath, outpath):
    # read the sample id to barcode correspondence
    lima_counts=pd.read_csv(lima_counts_by_moviename_path, sep='\t')
    lima_counts.sort_values(['kinnex_adapter', 'isoseq_primer'], ascending=[True, True], inplace=True)
    bam_path=os.path.join(bampath,"*.refine.bam")
    file_paths = sorted(glob.glob(bam_path))
    print(file_paths)
    refine_df=pd.DataFrame(columns=['movie_name','kinnex_adapter','isoseq_primer','refineBam'])
    for f in file_paths: 
        hifi_id=re.search(r'([\w.-]+).refine.bam', os.path.basename(f)).group(1)
        refine_list=[hifi_id.split('.')[0], hifi_id.split('.')[1], hifi_id.split('.')[2], f]
        refine_df.loc[len(refine_df)]=refine_list
    lima_counts=lima_counts.merge(refine_df, how='left', on=['movie_name','kinnex_adapter','isoseq_primer'])
    
    lima_counts[['sample_id', 'refineBam']].groupby('sample_id').apply(lambda bams: \
         merge_bams(bams.refineBam,str(bams.sample_id.unique()[0]), outpath) \
             if (np.all(pd.notnull(bams.refineBam))) else bams)
    return


def test_merge_bams(grouped_bams_list, op_bam_id, outpath): 
    print(grouped_bams_list)
    if os.path.exists(outpath):
        opath=os.path.join(outpath,str(op_bam_id+".merged.bam"))
    print(opath)
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-idmap",required=True, help="path to id map file")
    parser.add_argument("-bampath", required=False, help="path to input bams to be merged")
    parser.add_argument("-limacountsdir",required=True, help="path to lima counts file")
    parser.add_argument("-outdir", required=True, help="path to store outs") # create a folder called merged here- to add merged sorted bams
    parser.add_argument("-mergeReplicates", required=False, action='store_true', \
                        help="set to True if technical replicates are to be merged as described by the idmap file. \
                            If set to False only plot read counts from lima")
    parser.add_argument("-setTitleSamplePlot",required=False,\
                        help=" set plot title for readcounts by sample merged across replicates" ,\
                        default="Readcounts by sample - replicates combined" , type=str)

    args = parser.parse_args()

    if os.path.isfile(args.idmap):
        idmap_cleaned_df=sanity_idmap(args.idmap)
        if os.path.exists(args.limacountsdir):
            print(args.limacountsdir)
            countspaths=lima_cnt_combine(args.limacountsdir, args.outdir, idmap_cleaned_df)
            print(countspaths)
            if args.setTitleSamplePlot:
                plot_counts(countspaths, args.outdir, args.setTitleSamplePlot)
            else:
                plot_counts(countspaths, args.outdir)    
        else:
            print(f'ERROR: The file or directory at {args.limacountsdir} does not exist. Provide directory with lima demux outs.')
        if args.mergeReplicates and args.bampath :
            countsPathbyMovie=os.path.join(args.outdir,"lima_counts_by_moviename.tsv")
            if os.path.exists(args.bampath):
                opath=os.path.join(args.outdir,"merge") # create the out folder here for non-wdl version
                print('Making dir to add merged reads...')
                print(opath)
                os.mkdir(opath)
                group_bams(countsPathbyMovie,args.bampath,opath)
            else:
                print(f'ERROR: The file or directory at {args.bampath} does not exist. Please provide path to directory with bams to be merged here.') 
        else:
            print("No merging requested")         
    else:
        print(f'ERROR: The file or directory at {args.idmap} does not exist. Please provide idmap file here.')    
        

if __name__ == "__main__":
    main()
