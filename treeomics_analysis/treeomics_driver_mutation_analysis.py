import pandas as pd, os, glob

def treeomics_driver_mutation_vaf(idir,odir,pid):
    '''
    Extract driver mutations from treeomics output data

    RUN: python(idir = '/Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/PID_4/treeomics/',
                odir = '/path/to/output/dir/',
                pid = 'PID_4')
    
    '''
    df = pd.read_csv(glob.glob(f'{idir}/{pid}/treeomics/mut_*_variants.csv')[0],
                     sep=",",header=0)#, encoding='iso-8859-1')#read mutation variants treeomics output
    #pd.set_option('display.max_columns', None)
    df = df[df.Driver==True][['GeneSymbol', 'Driver', 'Phylogeny', 'BranchName'] + [col for col in df.columns if col.startswith('VAF')]]
    df.to_csv(f'{odir}/{pid}/{pid}.driver_mutation.vad.txt',sep="\t",header=True, index=False)
    print(f'Output generated:\n\t{odir}/{pid}/{pid}.driver_mutation.vad.txt')
    
    #return(df)


treeomics_driver_mutation_vaf(idir = '/Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/',
                              odir = '/Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/',
                              pid='PID_4')
