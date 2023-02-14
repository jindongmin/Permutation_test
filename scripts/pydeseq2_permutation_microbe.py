import pydeseq2

def pydeseq2_permutation_test(counts, group_labels, num_permutations=10):
    """Permutation test for differential expression using PyDESeq2.

    Parameters
    ----------
    counts : microbial count matrix, (n_samples, m_microbes).
    group_labels : array_like
        The group labels of shape (n_samples,).
    num_permutations : int, optional (default=1000)
        The number of permutations to perform.

    Returns
    -------
    p_values : array_like
        The p-values of the permutation test for each microbe.
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
    
    # Initialize the null distribution of the test statistic
    null_distribution = np.zeros((num_permutations, counts.shape[1]))

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

        # Add the permuted test statistic to the null distribution
        null_distribution[i, :] = permuted_stat_res.results_df['log2FoldChange']
    
    #transform the null distribution and real results to fold changes
    null_fold_changes = 2 ** null_distribution - 1
    real_fold_changes = 2 ** stat_res.results_df['log2FoldChange'].values - 1
    
    # Calculate the p-values for each gene
    p_values = np.zeros(counts.shape[1])
    for i in range(counts.shape[1]):
        if real_fold_changes[i] > 0:
            p_values[i] = np.mean(null_fold_changes[:, i] > real_fold_changes[i])
        else:
            p_values[i] = np.mean(null_fold_changes[:, i] < real_fold_changes[i])

    return p_values

    
    
#     # Calculate the p-values for each gene
#     p_values = np.zeros(counts.shape[1])
#     for i in range(counts.shape[1]):
#         if stat_res.results_df['log2FoldChange'].values[i] > 0:
#             p_values[i] = np.mean(null_distribution[:, i] > stat_res.results_df['log2FoldChange'].values[i])
#         else:
#             p_values[i] = np.mean(null_distribution[:, i] < stat_res.results_df['log2FoldChange'].values[i])

#     return p_values

