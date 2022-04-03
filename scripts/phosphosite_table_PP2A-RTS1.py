import os
import sys
import pandas as pd
import itertools


# Import phospho site tables
root = '../'
table_dir = os.path.join(root, 'tables')
setup_names = ['rts1D', 'SR', 'SR-hog1as', 'SR-rck2D']

tables = []
for setup_name in setup_names:
    path = os.path.join(table_dir, 'phospho_' + setup_name + '.tsv')
    tables.append(pd.read_csv(path, sep='\t', index_col=False))
combined_table = pd.concat(tables)


# Generate a pivote table of phosphorylation sites with mean, count, std 
pivoted_table = pd.pivot_table(
    data=combined_table, columns='Setup', values='Ratio Log2 normalized',
    aggfunc={'Ratio Log2 normalized': ['mean', 'count', 'std']},
    index=[
        'Protein Standard Name', 'Protein Systematic Name',
        'Phospho sites', 'Phospho (STY)',
    ],
)
pivoted_table.columns = pivoted_table.columns.swaplevel(0, 1)


# Replace missing and empty values
mean_fillna = '-'
count_fillna = 0
std_fillna = '-'
for setup in setup_names:
    pivoted_table[(setup, 'mean')].fillna(mean_fillna, inplace=True)
    pivoted_table[(setup, 'count')].fillna(count_fillna, inplace=True)
    pivoted_table[(setup, 'std')].fillna(std_fillna, inplace=True)


# Reorganize column order
col_order = itertools.chain(
    *[[(s, 'mean'), (s, 'count'), (s, 'std')] for s in setup_names]
)
pivoted_table = pivoted_table.ix[:, col_order]
pivoted_table.rename(
    columns={'mean': 'SILAC Ratio', 'count': 'Num', 'std': 'Std'}, inplace=True
)
pivoted_table.reset_index(inplace=True)


# Export table
path = os.path.join(root, 'tables', 'Supplementary_table_2_RTS1.tsv')
pivoted_table.to_csv(path, sep='\t', index=False)
