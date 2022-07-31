import os
import pandas as pd
def gencode_feature_coordinates(gtf, idir='./', odir='./', fout='out'):
    '''
    Extract FEATURE Specific Coordinates from Gencode GTF file (Region: PRI)
    
    RUN: gencode_feature_coordinates(idir='./', odir='./', fout = 'out.txt', gtf='gencode.gtf', feature) 

    Input:
            idir   : /path/to/gtf/
            odir   : /path/to/output/dir/
            fout   : out
            gtf    :  GTF file downloaded from Gencode (https://www.gencodegenes.org/human/)
                      GTF Format: https://www.gencodegenes.org/pages/data_format.html
            
            feature:  {'transcript' | 'exon' | 'gene' | 'start_codon' | 'stop_codon' | 'UTR' | 'CDS' | 'Selenocysteine'}

    
    Output: Returns Feature Coordinates and Gene Name Annotation
            Output file: 'out.feature.txt'
            chr     start   end     strand  symbol                                                info
            chr1    11869   12227   +       DDX11L1
            chr1    12613	12721   +       DDX11L1
    '''
    feature =  {1:'transcript', 2:'exon', 3:'gene', 4:'start_codon', 5:'stop_codon', 6:'UTR', 7:'CDS', 8:'Selenocysteine'}
    for k, v in feature.items(): print(f'"{k}": {v.upper()}')#print all the options 
    option = int(input('Type numeric option: '))
    if option in feature.keys():
        print('\nSelected feature:', feature[option])
        f = pd.read_csv(f'{idir}{gtf}', sep="\t", header=None, comment='#', usecols=[0,2,3,4,6,8], names=['chr','feature','start','end', 'strand', 'info'])#read specific columns from gtf file
        feat = f[f.feature == feature[option]]#[['chr','genomic_start','genomic_end', 'strand', 'info']]
        feat['symbol'] = feat['info'].str.split('gene_name "').str[1].str.split('";').str[0]
        feat[['chr','start','end', 'strand', 'symbol']].to_csv(f'{odir}{fout}.{option}.txt')
    else:
        print('Please enter correct Option'); exit()
    print(f'Successfully generated:\n\t{odir}{fout}.{option}.txt')
