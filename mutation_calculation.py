import os
import pandas as pd

from misc.sample_vs_gene_mut_table import funcotator_anot_vcf_mut_gene_table  as tb 

def main():
    '''
    Calculate Sample vs Gene Mutation Table
    '''
    path = '/Volumes/lilac_data_ziv/transciptome/paired_pnet/Project_12502_F/anotvcf/'
    sample = [f for f in os.listdir(path) if f.endswith('.pass.anot.vcf')]
    i = 0
    d = pd.DataFrame()
    for s in sample:
        if(i == 0):
            d = tb(f'{path}{s}')
            i += 1
        else:
            d = pd.concat([d,tb(f'{path}{s}')], axis=0)
        print(f'{s} is processed...')
    d.to_csv('sample_vs_gene_mut.table.txt', sep="\t", header=True, index=False)

        

if __name__ == '__main__':
    main()