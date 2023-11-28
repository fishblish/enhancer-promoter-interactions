# %%
import os
path = os.getcwd()
path = path.split('interactions', 1)[0] + 'interactions'
print(path)

# %%
import pandas as pd

bed = pd.read_csv(path+'/results/whole_data/out_params_all_chrom.bedpe', delimiter='\t', header=None)
bed.columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'x', 'score', 'y', 'z', 'color']
print(bed.shape)
print(bed['chr1'].unique())
bed.head(10)


# %%
bed['x'].replace('.', 0, inplace=True)
bed['y'].replace('.', 0, inplace=True)
bed['z'].replace('.', 0, inplace=True)
print(sum(bed['x']), sum(bed['y']), sum(bed['z']))
bed.drop(['x', 'y', 'z'], axis=1, inplace=True)
cicero_conns = bed.copy()
bed.head()

# %%
sum(cicero_conns['chr1'] != cicero_conns['chr2'])
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
chromatin_loops.shape

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
cicero_conns['intersection'] = ['n']*cicero_conns.shape[0]
cicero_conns.reset_index(level=0, inplace=True)
cicero_conns.head() 

# %%
cicero_conns = cicero_conns[cicero_conns['start1'] == 17260794]
results = {'chr': [], 'intersection1': [], 
                                'intersection2': [],
                                'cicero1': [], 'cicero2': [], 
                                'loop1': [], 'loop2': []}
for chrom in ['2L']:
    cicero = cicero_conns[cicero_conns['chr'] == chrom]
    loops = loop_conns[loop_conns['chr'] == chrom]
    for cicero_index, cicero_row in cicero.iterrows():
        for loop_index, loop_row in loops.iterrows():
            cicero1 = [cicero_row['start1'], cicero_row['end1']]
            cicero2 = [cicero_row['start2'], cicero_row['end2']]
            loops1 = [loop_row['range_start1'], loop_row['range_end1']]
            loops2 = [loop_row['range_start2'], loop_row['range_end2']]
            print(cicero1, cicero2, loops1, loops2)
            if (intersect_intervals(cicero1, loops1) and intersect_intervals(cicero2, loops2)):
                print('znalezione!')
                cicero_conns.loc[cicero_conns['index'] == cicero_row['index'], 'intersection'] = 'y'
                results['chr'].append(chrom)
                results['intersection1'].append(intersect_intervals(cicero1, loops1))
                results['intersection2'].append(intersect_intervals(cicero2, loops2))
                results['cicero1'].append(cicero1)
                results['cicero2'].append(cicero2)
                results['loop1'].append(loops1)
                results['loop2'].append(loops2)
                print(results)
            else:
                if intersect_intervals(cicero1, loops2) and intersect_intervals(cicero2, loops1):
                    print('Cos znalazalo!')
                    cicero_conns.loc[cicero_conns['index'] == cicero_row['index'], 'intersection'] = 'y'
                    results['chr'].append(chrom)
                    results['intersection1'].append(intersect_intervals(cicero1, loops2))
                    results['intersection2'].append(intersect_intervals(cicero2, loops1))
                    results['cicero1'].append(cicero1)
                    results['cicero2'].append(cicero2)
                    results['loop1'].append(loops1)
                    results['loop2'].append(loops2)
                    print(results)



# %%
