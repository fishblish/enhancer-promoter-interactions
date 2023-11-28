# %%
import os
path = os.getcwd()
path = path.split('interactions', 1)[0] + 'interactions'
print(path)

# %%
import pyreadr

result = pyreadr.read_r(path+'/results/whole_data/saved_conns_params.rds')
cicero_conns =result[None]
cicero_conns.head()

# %%
cicero_conns.shape

# %%
for i in ['1','2']:
    cicero_conns[['chr'+i, 'range'+i]] = cicero_conns['Peak'+i].str.split(':', 1, expand=True)
    cicero_conns[['start'+i, 'end'+i]] = cicero_conns['range'+i].str.split('-', 1, expand=True)
    cicero_conns = cicero_conns.drop(['Peak'+i, 'range'+i], axis=1)

cicero_conns.head()

# %%
sum(cicero_conns['chr1'] != cicero_conns['chr2'])

# %%
cicero_conns = cicero_conns.drop('chr2', axis=1)
cicero_conns = cicero_conns.rename(columns={'chr1': 'chr'})
cicero_conns.head()

# %%
import pandas as pd

chromatin_loops = pd.read_excel(path+'/data/muszka/TableS2_Dmel_loops.xlsx')
chromatin_loops.head()

# %%
unique_col = ['loop ID', 'loop size (in bp)', 'loop clustering', 'loop type']
columns_to_duplicate = [col for col in list(chromatin_loops.columns) if col not in unique_col+['loop ID']]
print(columns_to_duplicate)
col_dict = {}
for col in columns_to_duplicate:
    col_dict[col] = ['first', 'last']
for col in unique_col:
    col_dict[col] = 'first'

# %%
# join pairs of rows in chromatin_loops  with the same loop ID

print(chromatin_loops.shape)
chromatin_loops = chromatin_loops.groupby('loop ID').agg(col_dict)
chromatin_loops.head()


# %%
chromatin_loops.columns = chromatin_loops.columns.map('_'.join)
chromatin_loops.head()

# %%
# cols anchor chr_first and anchor chr_last are identical
sum(chromatin_loops['anchor chr_first'] != chromatin_loops['anchor chr_last'])

# %%
chromatin_loops.drop(columns=['anchor chr_last'], inplace=True)
chromatin_loops.rename(columns={'anchor chr_first':'chr', 
                                'anchor start_first':'start1', 
                                'anchor start_last':'start2',
                                'anchor end_first':'end1',
                                'anchor end_last':'end2'}, inplace=True)
loop_conns = chromatin_loops[['chr', 'start1', 'end1', 'start2', 'end2']].copy()
loop_conns.head()

# %%
for col in ['start1', 'start2']:
    loop_conns.loc[:,'range_'+ col] = loop_conns.loc[:,col]+5000

for col in ['end1', 'end2']:
    loop_conns.loc[:,'range_' + col] = loop_conns.loc[:,col].apply(lambda x: max(0, x - 5000))
loop_conns['chr'] = loop_conns['chr'].apply(lambda x: x[3:])

loop_conns.head()

# %%
cicero_conns.head()
#change type to int
cicero_conns[['start1', 'start2', 'end1', 'end2']] = cicero_conns[['start1', 'start2', 'end1', 'end2']].astype(int)

# %%
def intersect_intervals(interval1, interval2):
    inter_start = max(interval1[0], interval2[0])
    inter_end = min(interval1[1], interval2[1])
    if inter_start <= inter_end:
        return [inter_start, inter_end]
    else:
        return False

# %%
results = []
for chrom in set(cicero_conns['chr']).intersection(set(loop_conns['chr'])):
    cicero = cicero_conns[cicero_conns['chr'] == chrom]
    loops = loop_conns[loop_conns['chr'] == chrom]
    for (cicero1, cicero2) in zip(cicero[['start1', 'end1']].values, cicero[['start2', 'end2']].values):
        for (loops1, loops2) in zip(loops[['start1', 'end1']].values, loops[['start2', 'end2']].values):
            if intersect_intervals(cicero1, loops1) and intersect_intervals(cicero2, loops2):
                results.append({'chr': chrom, 'intersection1': intersect_intervals(cicero1, loops1), 
                                'intersection2': intersect_intervals(cicero2, loops2),
                                'cicero1': cicero1, 'cicero2': cicero2, 
                                'loop1': loops1, 'loop2': loops2})
                print(results[-1])
            else:
                if intersect_intervals(cicero1, loops2) and intersect_intervals(cicero2, loops1):
                    results.append({'chr': chrom, 'intersection1': intersect_intervals(cicero1, loops2), 
                                    'intersection2': intersect_intervals(cicero2, loops1),
                                    'cicero1': cicero1, 'cicero2': cicero2, 
                                    'loop1': loops2, 'loop2': loops1})
                    print(results[-1])





