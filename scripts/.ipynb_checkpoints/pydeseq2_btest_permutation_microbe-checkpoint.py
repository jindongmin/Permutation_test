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


def get_kegg_from_annotation(ids):
    """
    ids: Species_rep of the top/bottom k microbes
    """
    df_list = []
    for i in ids:
        f_name = "../../permutation_table/AD_species/{}_eggNOG.tsv".format(i)
        df_parsed = parse_genome(pd.read_table(f_name))
        df_list.append(df_parsed)
    df_cat = pd.concat(df_list, axis = 0)
    kegg_counts = to_sparse_matrix(df_cat)
    return kegg_counts
 
    
def _test(X, y, z, k, M, G, nc, pt):
    """
    X: count matrix of microbes (m microbes X n samples)
    y: disease lables for samples (n samples)
    z: colname of the condition (eg. disease status)
    k: top k microbes, eg k = 100
    M: metadata = pd.read_table('../table/eggNOG_species_rep.txt')
    G: annotation matrix of microbes (microbe ids X kegg ids)
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
    top_rep = M[M['Species'].isin(top_microbe)]
    bot_rep = M[M['Species'].isin(bot_microbe)]
    top_id = top_rep['Species_rep'].drop_duplicates()
    top_id = [_test_top_id for _test_top_id in top_id if "." not in _test_top_id]
    bot_id = bot_rep['Species_rep'].drop_duplicates()
    bot_id = [_test_bot_id for _test_bot_id in bot_id if "." not in _test_bot_id]
    #get top and bot G gene matrix
    top_gene_matrix = G.filter(items = top_id, axis=0)
    bot_gene_matrix = G.filter(items = bot_id, axis=0)
    frames = [top_gene_matrix, bot_gene_matrix]
    Conta_gene_matrix = pd.concat(frames)
    Conta_gene_matrix = Conta_gene_matrix.loc[:, (Conta_gene_matrix != 0).any(axis=0)]
    top_gene_matrix = Conta_gene_matrix.filter(items = top_id, axis=0)
    bot_gene_matrix = Conta_gene_matrix.filter(items = bot_id, axis=0)
    #binomial test Return number of genes elevated in top
    """
    C = btest(G[top], G[bot])
    return C
    """
    kegg = btest(top_gene_matrix,bot_gene_matrix,return_proportions=True)
    kegg_top = kegg.loc[kegg['side'] == 'groupA']
    kegg_top_sig = kegg_top.loc[kegg_top['pval'] <= pt]
    C = len(kegg_top_sig.index)
    return C    


def permutation_test(X, y, z, k, p, M, G, nc, pt):
    """
    p: number of permutations, eg p = 1000
    
    """
    T = _test(X, y, z, k, M, G, nc, pt)
    T_list = [0] * p
    for i in range(p):
        #shuffle the group lables
        y_permutated = np.random.permutation(y[z])
        y_permutated = pd.DataFrame(y_permutated, index=y.index)
        y_permutated.columns = [z]
        y_permutated.reindex(X.index)
        T_ = _test(X, y_permutated, z, k, M, G, nc, pt)
        T_list[i] = T_
        p_value = np.sum(T_list[i] > T) / (p+1)
    return T, p_value


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
    parser.add_argument("-M", dest="metadata",
                       help="metadata for ids and species",
                       required=True)
    parser.add_argument("-G", dest="gene_matrix",
                       help="gene matrix for microbes",
                       required=True)
    parser.add_argument("-nc", dest="number_cpus",default=8, type=int,
                       help="number of cpus")
    parser.add_argument("-pt", dest="pval",default=0.001, type=float,
                       help="threshold of p value")
    parser.add_argument("-o", dest="output_file",
                       help="number of genes and p value",
                       required=True)
    options = parser.parse_args()
    return options


def main():
    option = get_options()
    
    # input_1 = argv[1]
    # input_2 = argv[2]
    # input_3 = int(argv[3])
    # input_4 = int(argv[4])
    # input_5 = argv[5]
    # output_1 = argv[6]
    # output_2 = argv[5]

    input_table = pd.read_table(option.microbe_table, sep = '\t', index_col = 0)
    input_table2 = pd.read_table(option.disease_labels, sep = '\t', index_col = 'featureid')
    input_table3 = pd.read_table(option.metadata, sep = '\t')
    input_table4 = pd.read_table(option.gene_matrix, sep='\t', index_col = 0 )

    a,b = permutation_test(X = input_table,
                           y = input_table2,
                           z = option.disease_status,
                           k = option.top_k,
                           p = option.permutations,
                           M = input_table3,
                           G = input_table4,
                           nc = option.number_cpus,
                           pt = option.pval)

    with open(option.output_file, "w") as output_h:
        output_h.write(str(a) + "\n")
        output_h.write(str(b) + "\n")


if __name__ == '__main__':
    main()
