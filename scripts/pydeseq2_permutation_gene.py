import pydeseq2

def pydeseq2_permutation_test(counts, group_labels, num_permutations=10):
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
    refit_cooks=True#,
    #n_cpus=8,
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
    ranks = np.argsort(sorted_log2_fold_changes)

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
        permuted_stat_res = DeseqStats(permuted_dds)#, n_cpus=8)
        # Get the test statistics for each gene with the permuted group labels
        permuted_res = permuted_stat_res.summary()
        log2_fold_changes_permutated = permuted_stat_res.results_df["log2FoldChange"]
        sorted_indices_permutated = np.argsort(log2_fold_changes_permutated)[::-1]
        sorted_log2_fold_changes_permutated = log2_fold_changes_permutated[sorted_indices_permutated]
        ranks_permutated_ori = np.argsort(sorted_log2_fold_changes_permutated)
        ranks_permutated = ranks_permutated_ori.reindex(ranks.index)
        #add ranks to the null distribution
        null_distribution[i,:] = ranks_permutated

    #TO DO: calculate the p values based on ranks
    p_values_rank = np.zeros(counts.shape[1])
    for i in range(counts.shape[1]):
        p_values_rank[i] = np.sum(null_distribution[:, i] > ranks[i]) / num_permutations

    return p_values_rank, stat_res.results_df
