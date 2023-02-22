#! /usr/bin/env python

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import numpy as np
import pydeseq2
from sys import argv
import pandas as pd

def pydeseq2_permutation_test(counts, group_labels, num_permutations=1000):
    """Permutation test for differential expression using PyDESeq2.

    Parameters
    ----------
    counts : array_like
        The count matrix of shape (n_samples, n_genes).
    group_labels : array_like
        The group labels of shape (n_samples,).
    num_permutations : int, optional (default=1000)
        The number of permutations to perform.

    Returns
    -------
    p_values : array_like
        The p-values of the permutation test for each gene.
    """
    # Convert counts and group labels to PyDESeq2 input format
    dds = DeseqDataSet(
    counts,
    group_labels,
    design_factors="disease",  # compare samples based on the "disease"
    # column ("B" vs "A")
    refit_cooks=True,
    n_cpus=8
)
    dds.deseq2()
    stat_res = DeseqStats(dds, n_cpus=8)
    # Get the test statistics for each gene
    res = stat_res.summary()
    res_df = stat_res.results_df
    sorted_df = res_df.sort_values(by='log2FoldChange', ascending=False)
    rankings = sorted_df['log2FoldChange'].rank(ascending=False)
    sorted_df['Rankings'] = rankings
    sorted_df.dropna()
    sorted_df = sorted_df.sort_index()
    ranks = sorted_df['Rankings']
    #initialize the null distribution of the ranking
    null_distribution = np.zeros((num_permutations, counts.shape[1]))

    # Perform the permutations
    for i in range(num_permutations):
        # Shuffle the group labels
        permuted_labels = np.random.permutation(group_labels['disease'])
        permuted_labels = pd.DataFrame(permuted_labels, index=group_labels.index)
        permuted_labels = pd.DataFrame(permuted_labels)
        permuted_labels.columns = ['disease']
        # Convert the permuted group labels to PyDESeq2 input format
        # Run the DESeq2 analysis with the permuted group labels
        permuted_labels.reindex(counts.index)
        #counts.reindex(permuted_labels.index)
        permuted_dds = DeseqDataSet(
        counts, 
        permuted_labels, 
        design_factors="disease")
        permuted_dds.deseq2()
        permuted_stat_res = DeseqStats(permuted_dds, n_cpus=8)
        # Get the test statistics for each gene with the permuted group labels
        permuted_res = permuted_stat_res.summary()
        res_df_permuted = permuted_stat_res.results_df
        sorted_df_permuted = res_df_permuted.sort_values(by='log2FoldChange', ascending=False)
        rankings_permuted = sorted_df_permuted['log2FoldChange'].rank(ascending=False)
        sorted_df_permuted['Rankings'] = rankings_permuted
        sorted_df_permuted.dropna()
        sorted_df_permuted = sorted_df_permuted.sort_index()
        #add ranks to the null distribution
        null_distribution[i,:] = sorted_df_permuted['Rankings'].values

    #TO DO: calculate the p values based on ranks
    p_values_rank = np.zeros(counts.shape[1])
    for i in range(counts.shape[1]):
        p_values_rank[i] = np.sum(null_distribution[:, i] > ranks[i]) / num_permutations

    return p_values_rank, sorted_df

input_1 = argv[1]
input_2 = argv[2]
output_1 = argv[3]
output_2 = argv[4]

input_table = pd.read_table(input_1, sep = '\t', index_col = 0)
input_table2 = pd.read_table(input_2, sep = '\t', index_col = 0)
a,b = pydeseq2_permutation_test(input_table,input_table2)
#a.to_csv(output_1)
np.savetxt(output_1,a,delimiter=",")
b.to_csv(output_2, sep = '\t')
