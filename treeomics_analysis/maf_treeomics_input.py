import pandas as pd

def convert_maf_treeomics_input(maf_file):
    '''
    Treeomics Input files:
        1. Coverage.txt
            Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
        2. mut_read_table.txt
            Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
    '''
    cols = ['Chromosome','Start_Position','Tumor_Seq_Allele1','Tumor_Seq_Allele2','Hugo_Symbol','Tumor_Sample_Barcode',
        't_depth','t_alt_count']#coverage.txt, mut_read_table.txt
    maf = pd.read_csv(maf_file, sep="\t", header=0, usecols=cols, comment='#')
    maf['Change'] = maf.Tumor_Seq_Allele1 + '>' + maf.Tumor_Seq_Allele2
    maf = maf.rename(columns={'Hugo_Symbol':'Gene', 'Start_Position':'Position'})
    
    #coverage
    cov = maf[['Chromosome','Position', 'Change', 'Gene', 'Tumor_Sample_Barcode','t_depth']]
    cov = cov.pivot_table(index=['Chromosome', 'Position', 'Change', 'Gene'],
                              columns='Tumor_Sample_Barcode',
                              values='t_depth',
                              aggfunc='first')
    cov = cov.reset_index()# Reset the index
    cov.columns.name = None# Rename the columns
    print(cov.head())
    cov.to_csv('coverage.txt',sep="\t",header=True, index=False)

    #mut_read_table
    mut = maf[['Chromosome','Position', 'Change', 'Gene', 'Tumor_Sample_Barcode','t_alt_count']]
    mut = mut.pivot_table(index=['Chromosome', 'Position', 'Change', 'Gene'],
                              columns='Tumor_Sample_Barcode',
                              values='t_alt_count',
                              aggfunc='first')
    mut = mut.reset_index()# Reset the index
    mut.columns.name = None# Rename the columns
    mut.to_csv('mut_read_table.txt',sep="\t",header=True, index=False)
    print(mut.head())

convert_maf_treeomics_input(maf_file='Proj_08138_C___FILLOUT.V3b.maf.txt')
