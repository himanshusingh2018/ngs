import pandas as pd
def funcotator_anot_vcf_mut_gene_table(vcf):
    '''
    RUN: funcotator_anot_vcf_mut_gene_table('file.vcf')
    Input: Funcotator VCF file
    Output: 
        fname   nGene   nMutation
        file1   14      20
        file2   18      30
    '''
    df = pd.read_csv(vcf, sep="\t", header=0, skiprows=266, usecols=['INFO'], encoding='latin-1')
    df = pd.DataFrame(df.INFO.str.split('\=\[').str[1].str.split('|').str[0])
    data = [vcf.split('/')[-1], df.size, len(df.drop_duplicates()), len(df[df.INFO=='Unknown'])]
    table = pd.DataFrame([data],
                         columns=['sample', 'num_genes', 'num_mutation', 'unannotated_mutation'])
    return(table)