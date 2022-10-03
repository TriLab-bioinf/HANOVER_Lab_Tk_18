#!/home/lorenziha/data/miniconda3/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import re
import collections
from numpy import size
import numpy as np
from tqdm import tqdm
import argparse
from statistics import median

######################################################################################
# Process parameters
######################################################################################

# Initialize parser
parser = argparse.ArgumentParser(
    description="This proram process a methylation bed file (zero-based coords) to generate methylation distribution plots across all genes in a genome\n\n"
)

# Input should look like:
# chr1	7157890	7157891	-	0	4	CHH	CAA	ENSMUSG00000051285
# Fields = chromosome, meth-start, meth-end, strand, # meth reads, # unmeth reads, methyl type (it uses CG only), context sequence, gene_id that matches gtf annotation file
#
# Command example:
# python ./plot_methyl_distrib_across_genes.py -g ./Mus_musculus.GRCm39.107.chr.gtf -m oga1.meth_geneid.bed -G OGA -o oga1

# Adding argument
parser.add_argument(
    "-m",
    "--methylation",
    help="methylation bed file",
    type=str,
    required=False,
)

parser.add_argument(
    "-g",
    "--annotation_gtf",
    help="annotation in gtf format (gene IDs must match those from methylation bed file)",
    type=str,
    required=False,
)

parser.add_argument(
    "-G",
    "--genotype",
    help="genotype  of cell line. Used to group samples and generate plots (e.g. OGA or WT).",
    type=str,
    required=False,
    default='NONE',
)

parser.add_argument(
    "-o",
    "--output_prefix",
    help="output prefix",
    type=str,
    required=False,
)

args = parser.parse_args()
methFile = args.methylation
annotFile = args.annotation_gtf
outputPrefix = args.output_prefix
my_genotype = args.genotype

######################################################################################
# Define functions
######################################################################################

def get_relative_coords(my_meth_df):
    flank = 2000
    resolution = 1 # number of digits in the normalized gene coords (0=lowest resolution)
    my_df = pd.DataFrame( columns=['chr','rel_pos','meth','flag','group']) #,'target','end5','end3','strand','geneid']) # Stores methylation sites with gene-relative position and % methylation per site per gene per sample

    # Iterate thought filtered methylation sites and convert coords into gene-relative coords (0% => TSS; 100% => TES)
    #for i in tqdm(my_meth_df.loc[0:4000].index):
    for i in tqdm(my_meth_df.index):
        
        # Extract info from methylation df
        my_group = my_meth_df.loc[i, "group"]
        my_geneid = my_meth_df.loc[i, "geneid"]
        my_target = my_meth_df.loc[i, "start"]
        my_perc = my_meth_df.loc[i, "meth"] * 100 /  (my_meth_df.loc[i, "not_meth"] + my_meth_df.loc[i, "meth"])
        my_chr = my_meth_df.loc[i, "chr"]

        # Extract gene data from gff dicts
        my_start = gff_end5[my_geneid]
        my_end = gff_end3[my_geneid]
        my_strand = gff_strand[my_geneid]

        # Coord -> Relative coord
        if (my_target >= my_start) & (my_target <= my_end):
            # CODING
            if my_strand == '+':
                rel_pos = round(100 * (my_target - my_start) / (my_end - my_start),resolution)
                my_flag = 'coding'
            else:
                rel_pos = round(100 * (my_end - my_target) / (my_end - my_start),resolution)
                my_flag = 'coding'

        elif (my_target >= my_start - flank) & (my_target <= my_start):
            # UPSTREAM
            if my_strand == '+':
                rel_pos = round(100 * (my_target - my_start) / flank / 3,resolution)
                my_flag = 'upstream'
            else:
                rel_pos = round(100 + (100 * (my_start - my_target) / flank / 3), resolution)
                my_flag = 'downstream'

        elif (my_target >= my_end) & (my_target <= my_end + flank):
            # DOWNSTREAM
            if my_strand == '+':
                rel_pos = round(100 + 100 * ((my_target - my_end) / flank / 3),resolution)
                my_flag = 'downstream'
            else:
                rel_pos = round(100 * (my_end - my_target) / flank / 3, resolution)
                my_flag = 'upstream'

        my_df = my_df.append({'chr':my_chr, 'rel_pos':rel_pos, 'meth':my_perc,'flag':my_flag,'group':my_group}, ignore_index=True)

    return(my_df)


######################################################################################
# Main
######################################################################################

# Read gff3 file
print(f'Loading annotation gtf file ... {annotFile}')
gff_df = pd.read_csv(annotFile, sep="\t", comment='#', names=('chr','flag', 'feat','end5','end3','score','strand','phase','desc'), dtype={'chr':str})

# Keep only gene info and protein_coding
gene_gff = gff_df[gff_df['feat'] == 'gene']
gene_gff = gene_gff.loc[gene_gff['desc'].str.contains( pat="protein_coding")]

# Adding geneid column
my_gene_id = []
for gene in gene_gff['desc']:
    gene_desc = re.search("(ENSMUSG\d+)", gene)
    my_gene_id.append(gene_desc[0])

gene_gff.loc[:,'geneid'] = my_gene_id

# Make dictionaries for gene coords and strand
gff_end5 = dict(zip(tuple(gene_gff['geneid']),gene_gff['end5']))
gff_end3 = dict(zip(tuple(gene_gff['geneid']),gene_gff['end3']))
gff_strand = dict(zip(tuple(gene_gff['geneid']),gene_gff['strand']))

# Load methylation data
print(f'Lodaing methylation file ... {methFile}')
my_meth_df = pd.read_csv(methFile, sep="\t", comment='#', header=None, names=['chr','start','end','strand','meth','not_meth','type','seq', 'geneid'] )

# Adding genotype column
my_meth_df["group"] = my_genotype

# Keeping just CG methypation sites
print("Keeping only CG type")
my_meth_df = my_meth_df.loc[my_meth_df["type"] == "CG"]

# Converting absolute pos to relative pos
print(f"Getting relative position for {methFile}")
my_rel_meth_df = get_relative_coords(my_meth_df=my_meth_df)
        
# write my_rel_meth_df to a file
my_csv_file = outputPrefix + '.relPos.csv'
print(f"Writing relative position file ... {my_csv_file}")
my_rel_meth_df.to_csv(my_csv_file, index=False)
