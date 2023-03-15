#! /usr/bin/env python

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import random
from statsmodels.stats.proportion import binom_test
from statsmodels.stats.multitest import multipletests
#from sys import argv
from argparse import ArgumentParser 
import os


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
    #print("min n", min(list(n.values)))
    obs = list(zip(list(pa1.values), list((pa2.values + 1) / (pa2 + 1).sum()), list(n.values)))
    pvals = pd.Series([binom_test(_a, _n, _b, 'two-sided') for (_a, _b, _n) in obs],
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


def to_sparse_matrix(func_df, genome_id='genome_id', kegg_id='KEGG_ko'):
    # create genome-specific index
    ogus = list(set(func_df[genome_id]))
    ogu_lookup = pd.Series(np.arange(0, len(ogus)), ogus)
    # create KEGG-specific index
    keggs = list(set(func_df[kegg_id]))
    kegg_lookup = pd.Series(np.arange(0, len(keggs)), keggs)
    # rename names as numbers
    ogu_id = func_df[genome_id].apply(lambda x: ogu_lookup.loc[x]).astype(np.int64)
    kegg_id = func_df[kegg_id].apply(lambda x: kegg_lookup.loc[x]).astype(np.int64)
    # assign the presence / absence of a gene
    func_df['count'] = 1
    c = func_df['count'].values
    # format into a matrix
    data = coo_matrix((c, (ogu_id, kegg_id)))
    ko_ogu = pd.DataFrame(data.todense(), index=ogus, columns=keggs)
    return ko_ogu
 
    
def _test(X, y, z, k, G, nc, pt):
    """
    X: count matrix of microbes (m microbes X n samples)
    y: disease lables for samples (n samples)
    z: colname of the condition (eg. disease status)
    k: top k microbes, eg k = 100
    G: annotation matrix of microbes (microbe ids X kegg ids), and have a column with species names
    nc: number of cpu
    pt: the threshold of p value, eg pt = 0.001
    """
    # Convert counts and group labels to PyDESeq2 input format
    dds = DeseqDataSet(
        X, 
        y, 
        design_factors = z,# compare samples based on the "disease"
        refit_cooks = True,
        n_cpus = nc
    )
    dds.deseq2()
    #set of the sample status: Healthy, disease
    status = set(y[z])
    #remove healthy
    status.remove("Healthy")
    disease = ', '.join(status)
    #sample status colname: disease, 
    #disease: name of that disease
    stat_res = DeseqStats(dds,contrast = [z, disease, "Healthy" ], n_cpus=nc)
    #stat_res = DeseqStats(dds,contrast = ["disease", "Healthy", "AD" ], n_cpus=8)
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
    top_gene_matrix = G[G['Species'].isin(top_microbe)]
    bot_gene_matrix = G[G['Species'].isin(bot_microbe)]
    frames = [top_gene_matrix, bot_gene_matrix]
    #concatenate two gene matrix and filter out the kegg with only 0 values
    Conta_gene_matrix = pd.concat(frames)
    Conta_gene_matrix = Conta_gene_matrix.loc[:, (Conta_gene_matrix != 0).any(axis=0)]
    top_gene_matrix = Conta_gene_matrix[Conta_gene_matrix['Species'].isin(top_microbe)]
    bot_gene_matrix = Conta_gene_matrix[Conta_gene_matrix['Species'].isin(bot_microbe)]
    top_gene_matrix = top_gene_matrix.drop('Species', axis=1)
    bot_gene_matrix = bot_gene_matrix.drop('Species', axis=1)
    #binomial test Return number of genes elevated in top: diseased
    #C = btest(G[top], G[bot])
    #kegg = btest(top_gene_matrix,bot_gene_matrix,return_proportions=True)
    kegg = btest(bot_gene_matrix,top_gene_matrix,return_proportions=True)
    kegg_top = kegg.loc[kegg['side'] == 'groupB']
    kegg_top_sig = kegg_top.loc[kegg_top['pval'] <= pt]
    C = len(kegg_top_sig.index)
    return C, kegg_top_sig   


def permutation_test(X, y, z, k, p, G, nc, pt):
    #p: number of permutations, eg p = 1000
    T, kegg_top_sig = _test(X, y, z, k, G, nc, pt)
    print(kegg_top_sig)
    T_list = np.zeros(p)
    for i in range(p):
        #shuffle the group lables
        y_permutated = np.random.permutation(y[z])
        y_permutated = pd.DataFrame(y_permutated, index=y.index)
        y_permutated.columns = [z]
        y_permutated.reindex(X.index)
        y = y_permutated
        T_, foo = _test(X, y, z, k, G, nc, pt)
        T_list[i] = T_
    p_value = np.sum(T_list > T) / (p+1)
    return T, p_value, kegg_top_sig


def get_options():
    parser = ArgumentParser()
    parser.add_argument("-x", dest="microbe_table",
                       help="count matrix of microbes (m microbes X n samples)",
                       required=True)
    parser.add_argument("-y", dest="disease_labels",
                       help="disease lables for samples (n samples)",
                       required=True)
    parser.add_argument("-z", dest="disease_status",
                       help="colnames of the diseases status",
                       required=True)
    parser.add_argument("-k", dest="top_k", default=100, type=int,
                       help="top k microbes. "
                            "Default: %(default)i")
    parser.add_argument("-p", dest="permutations", default=1000, type=int,
                       help="number of permutations. "
                            "Default: %(default)i")
    parser.add_argument("-G", dest="gene_matrix",
                       help="gene matrix for microbes",
                       required=True)
    parser.add_argument("-nc", dest="number_cpus",default=8, type=int,
                       help="number of cpus")
    parser.add_argument("-pt", dest="pval",default=0.05, type=float,
                       help="threshold of p value")
    parser.add_argument("-o", dest="output_dir",
                       help="new output folder containing following files: "
                            "1) txt: number of genes and p value, "
                            "2) table: kegg_top_sig ",
                       required=True)
    options = parser.parse_args()
    os.mkdir(options.output_dir)
    return options


def main():
    option = get_options()
    input_table = pd.read_table(option.microbe_table, sep = '\t', index_col = 0)
    input_table2 = pd.read_table(option.disease_labels, sep = '\t', index_col = 'featureid')
    input_table3 = pd.read_table(option.gene_matrix, sep='\t', index_col = 0 )

    T, p_value, kegg_top_sig = permutation_test(X = input_table,
                                                y = input_table2,
                                                z = option.disease_status,
                                                k = option.top_k,
                                                p = option.permutations,
                                                G = input_table3,
                                                nc = option.number_cpus,
                                                pt = option.pval)

    with open(option.output_dir + "/n_genes_p_value.txt", "w") as output_h:
        output_h.write(str(T) + "\n")
        output_h.write(str(p_value) + "\n")
    kegg_top_sig.to_csv(option.output_dir + "/" + "kegg_sig.tsv", sep = '\t')


if __name__ == '__main__':
    main()
