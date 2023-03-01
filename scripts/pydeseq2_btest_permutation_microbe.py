#! /usr/bin/env python

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import random
from statsmodels.stats.proportion import binom_test
from statsmodels.stats.multitest import multipletests

def btest(pa1, pa2, seed=0, return_proportions=False):
    """ Performs genome wide binomial test between two groups of taxa
    Parameters
    ----------
    df1 : pd.DataFrame
        Rows are taxa, columns are genes
    df2 : pd.DataFrame
        Rows are taxa, columns are genes
    Returns
    -------
    pd.Series : list of genes associated with df1
    pd.Series : list of genes associated with df2
    """
    np.random.seed(seed)
    random.seed(seed)
    #pa1 = df1 > 0
    #pa2 = df2 > 0
    idx = list(set(pa1.columns) | set(pa2.columns))
    idx.sort()
    pa1 = pa1.sum(axis=0).reindex(idx).fillna(0)
    pa2 = pa2.sum(axis=0).reindex(idx).fillna(0)
    n = pa1 + pa2
    obs = list(zip(list(pa1.values), list((pa2.values + 1) / (pa2 + 1).sum()), list(n.values)))
    pvals = pd.Series([binom_test(a, n, b, 'two-sided') for (a, b, n) in obs],
                      index=n.index)
    if return_proportions:
        res = pd.DataFrame({'groupA': pa1, 'groupB': pa2, 'pval': pvals})
        def relabel_f(x):
            if x['groupA'] < x['groupB']:
                return 'groupB'
            else:
                return 'groupA'
        res['side'] = res.apply(relabel_f, axis=1)
        return res

    return pvals

def parse_genome(df):
    genome_id = df['#query'][0].split('_')[0]
    keggs = df['KEGG_ko'].replace('-', None).dropna()
    keggs = list(map(lambda x: x.split(','), keggs.values))
    keggs = sum(keggs, [])
    keggs = pd.DataFrame({'KEGG_ko': keggs})
    keggs['genome_id'] = genome_id
    return keggs

def get_kegg_from_annotation(ids):
    """
    ids: Species_rep of the top/bottom k microbes
    """
    df_list = []
    for i in ids:
        f_name = "../table/Species/{}_eggNOG.tsv".format(i)
        df_parsed = parse_genome(pd.read_table(f_name))
        df_list.append(df_parsed)
    df_cat = pd.concat(df_list, axis = 0)
    kegg_counts = to_sparse_matrix(df_cat)
    return kegg_counts
    
def _test(X, y, k, p, M):
    """
    X: count matrix of microbes (m microbes X n samples)
    y: disease lables for samples (n samples)
    G: gene table
    k: top k microbes, eg k = 100
    p: number of permutations, eg p = 1000
    M: metadata = pd.read_table('../table/eggNOG_species_rep.txt')
    """
    # Convert counts and group labels to PyDESeq2 input format
    dds = DeseqDataSet(
        X, 
        y, 
        design_factors = "disease", # compare samples based on the "disease"
        refit_cooks = True,
        n_cpus = 8
    )
    dds.deseq2()
    stat_res = DeseqStats(dds, n_cpus=8)
    res = stat_res.summary()
    res_df = stat_res.results_df
    res_df["CI_95"] = res_df["log2FoldChange"] + res_df['lfcSE'] * 1.96
    res_df["CI_5" ] = res_df["log2FoldChange"] - res_df['lfcSE'] * 1.96
    #top k microbes: the ones with largest CI_5
    top = res_df.sort_values(by=['CI_5'],ascending=False).head(k)
    #bottom k microbes: the ones with smallest CI_95
    bot = res_df.sort_values(by=['CI_95'],ascending=True).head(k)
    #convert microbe names to species ids
    top_microbe = set(top.index)
    bot_microbe = set(bot.index)
    top_rep = M[M['Species'].isin(top_microbe)]
    bot_rep = M[M['Species'].isin(bot_microbe)]
    top_id = top_rep['Species_rep'].drop_duplicates()
    bot_id = bot_rep['Species_rep'].drop_duplicates()
    #wget all ids: dir: /mnt/home/djin/ceph/Permutation_test/table/Species
    #get G gene matrix
    top_gene_matrix = get_kegg_from_annotation(top_id)
    bot_gene_matrix = get_kegg_from_annotation(bot_id)

    #binomial test Return number of genes elevated in top
    """
    C = btest(G[top], G[bot])
    return C
    """
    kegg = btest(top_gene_matrix,bot_gene_matrix,return_proportions=True)
    kegg_top = kegg.loc[kegg['side'] == 'groupA']
    kegg_top_sig = kegg_top.loc[kegg_top['pval'] <= 0.001]
    C = len(kegg_top_sig.index)
    return C    

def permutation_test(X, y, k, p, M):
    T = _test(X, y, k, p, M)
    for i in range(p):
        #shuffle the group lables
        y_permutated = np.random.permutation(y['disease'])
        y_permutated = pd.DataFrame(y_permutated, index=y.index)
        y_permutated.columns = ['disease']
        y_permutated.reindex(counts.index)
        T_ = _test(X, y_permutated, k, p, M)
    return T, (T > T_) / (p)
    #return T, (T > T_) / (p+1)





