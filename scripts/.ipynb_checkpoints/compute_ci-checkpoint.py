

def compute_ci(table_lfc):
    table_lfc['CI_5'] = table_lfc['log2FoldChange'] - table_lfc['lfcSE']*1.96
    table_lfc['CI_95'] = table_lfc['log2FoldChange'] + table_lfc['lfcSE']*1.96
    i_negative = table_lfc.sort_values(by=['CI_95'],ascending=True).head(2)
    i_positive = table_lfc.sort_values(by=['CI_5'],ascending=False).head(2)
    return i_negative, i_positive
