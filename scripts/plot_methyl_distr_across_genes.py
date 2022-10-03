#!/home/lorenziha/data/miniconda3/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import re
import collections
from numpy import size
import numpy as np
import argparse

######################################################################################
# Process parameters
######################################################################################

# Initialize parser
parser = argparse.ArgumentParser(
    description="This script process relPos methylation csv file to generate methylation distribution plots across all genes in a genome\n\n"
)

# Input should look like:
# chr,rel_pos,meth,flag,group
# chr1,-15.8,75.0,upstream,OGA
# chr1,-15.7,50.0,upstream,OGA

# Command example:
# 

# Adding argument
parser.add_argument(
    "-L",
    "--file_list",
    help="Comma-separated list of methylation relPos input files",
    type=str,
    required=False,
)

parser.add_argument(
    "-o",
    "--output_prefix",
    help="output prefix for plot files",
    type=str,
    required=False,
)

args = parser.parse_args()
methFiles = args.file_list
outputPrefix = args.output_prefix

######################################################################################
# Define functions
######################################################################################

def slide_window(_rel_pos, _methyl, k):
	my_res = pd.DataFrame(columns=['relPos','medMethyl'])
	_rel_pos = [int(x) for x in _rel_pos]
	mydict = dict(zip(_rel_pos, _methyl))
	my_gene_coords = list(np.arange(min(_rel_pos),max(_rel_pos),1))
	#print(my_gene_coords)
	n = len(my_gene_coords)

	for i in range(n-k):
		mykeys = my_gene_coords[i:i+k]
		#print(mykeys)
		my_wind_pos = np.mean(mykeys)
		#print(mykeys, my_wind_pos)
		# Get methylation values from keys
		myvals = [mydict.get(x,None) for x in mykeys]
		# Remove None from the window before estimating the median
		myvals = [i for i in myvals if i is not None]
		if len(myvals) == 0:
			continue
		
		mymed = np.mean(myvals)
		my_res = my_res.append(dict({'relPos':my_wind_pos, 'medMethyl':mymed}), ignore_index=True)
		
	return my_res

def calculate_mean_perc(my_meth_df):
    # Calculate mean perc methylation per methylation site across all genes and samples, within groups (OGA, WT)
    my_mean_meth = my_meth_df.groupby(['group','rel_pos'], as_index=False)['meth'].mean()
    return(my_mean_meth)

def plot_methylation_per_site(my_mean_meth):
    # Generate methylation distribution plot across genes
    groups = my_mean_meth['group'].unique()
    fig, ax = plt.subplots(figsize=(24, 5))
    
    # We need to generate one plot per group (OGA, WT) to assign different line colors
    for rank, my_group in enumerate(groups):
        my_x_vals = my_mean_meth.loc[my_mean_meth.group == my_group, 'rel_pos']
        my_y_vals = my_mean_meth.loc[my_mean_meth.group == my_group, 'meth']
        line = plt.plot(
            my_x_vals,
            my_y_vals,
            color=['red','blue'][rank], 
            alpha = 0.5, 
            label=my_group)
            
    plt.legend()
    #plt.show()
    plt.savefig('methyl_distribution_per_pos_genes.png')

def plot_methylation_per_window(my_results):
    # Generate methylation distribution plot across genes
    groups = my_results['group'].unique()
    fig, ax = plt.subplots(figsize=(24, 5))
    # We need to generate one plot per group (OGA, WT) to assign different line colors
    for rank, my_group in enumerate(groups):
        my_x_vals = my_results.loc[my_results.group == my_group, 'relPos']
        my_y_vals = my_results.loc[my_results.group == my_group, 'medMethyl']
        line = plt.plot(
            my_x_vals,
            my_y_vals,
            color=['red','blue'][rank], 
            alpha = 0.5, 
            label=my_group)
            
    plt.legend()
    #plt.show()
    plt.savefig('methyl_distribution_per_window_genes.png')

######################################################################################
# Main
######################################################################################

# Load methylation input files and concatenate them in a single pandas' data frame
my_file_list = list(methFiles.split(','))

meth_rel_all_df = pd.DataFrame( columns=['chr','rel_pos','meth','flag','group'])

for my_csv_file in my_file_list:
    my_rel_meth_df = pd.read_csv(my_csv_file, comment='#', header=1, names=['chr','rel_pos','meth','flag','group'])
    meth_rel_all_df = meth_rel_all_df.append(my_rel_meth_df)

# Sort mthyl sites by relPos
meth_rel_all_df = meth_rel_all_df.sort_values('rel_pos')

# Calculate mean perc methylation per methylation site across all genes and samples, within groups (OGA, WT)
mean_meth_df = calculate_mean_perc(meth_rel_all_df)
print('Generating per site methylation distribution plot')
plot_methylation_per_site(mean_meth_df)

# Calculate mean perc methylation per methylation site across all genes and samples, within groups (OGA, WT)
window_size = 10
window_results = pd.DataFrame(columns=['relPos','medMethyl','group'])
my_groups = list(meth_rel_all_df['group'].unique())

print('Calculating mean methylation per window')
for my_group in my_groups:
    my_group_met_df = meth_rel_all_df.loc[meth_rel_all_df['group'] == my_group]
    my_group_met_res = slide_window(_rel_pos=my_group_met_df['rel_pos'], _methyl=my_group_met_df['meth'], k=window_size)
    my_group_met_res['group'] = my_group
    window_results = window_results.append(my_group_met_res)

print('Generating per window methylation distribution plot')
plot_methylation_per_window(window_results)


