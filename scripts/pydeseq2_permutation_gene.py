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
    """Permutation test for differential abundance using PyDESeq2.

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
    n_cpus=8,
)
    dds.deseq2()
    stat_res = DeseqStats(dds, n_cpus=8)
    # Get the test statistics for each gene
    res = stat_res.summary()
    
    #extract the log2FoldChange for each feature
    log2_fold_changes = stat_res.results_df["log2FoldChange"]
    #sort the log2foldchanges in descending order
    sorted_indices = np.argsort(log2_fold_changes)[::-1]
    sorted_log2_fold_changes = log2_fold_changes[sorted_indices]
    #calculate the rank of each gene
    ranks = np.argsort(sorted_indices)


    # Perform the permutations
    for i in range(num_permutations):
        # Shuffle the group labels
#         permuted_labels = np.random.permutation(group_labels)
#         permuted_labels = pd.DataFrame(permuted_labels)  
        permuted_labels = np.random.permutation(group_labels['disease'])
        permuted_labels = pd.DataFrame(permuted_labels, index=group_labels.index)
        permuted_labels = pd.DataFrame(permuted_labels)
        permuted_labels.columns = ['disease']
        # Convert the permuted group labels to PyDESeq2 input format
        # Run the DESeq2 analysis with the permuted group labels
        counts.reindex(permuted_labels.index)
        permuted_dds = DeseqDataSet(
        counts, 
        permuted_labels, 
        design_factors="disease")
        permuted_dds.deseq2()
        permuted_stat_res = DeseqStats(permuted_dds, n_cpus=8)
        # Get the test statistics for each gene with the permuted group labels
        permuted_res = permuted_stat_res.summary()


    #p values based on ranks
    return p_values_ranks, stat_res.results_df

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
