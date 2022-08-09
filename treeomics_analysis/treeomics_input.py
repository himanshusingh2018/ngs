import os
import pandas as pd

def treeomics_input_vcfanotfuncotator(vcfdir='vcf', odir='./'):
    '''
    It generate input files from the VCF files from multiple samples located in one directory. 
    
    COMMAND:
    treeomics_input(vcfdir = './vcf/, odir='./)

    INPUT: Directory with all Funcotator Annotated GZ VCF Files

    OUTPUT:
    Generate Two Input files for TREEOMICS from Annotated VCF file

    1. <mut_reads.txt> mut-reads: number of reads reporting a variant (row) in each sample
       AD (Allele Depth) in VCF file (e.g. AD: 1,15; 1: ref and 15: alt)
       AD = 15

            mut-reads table: tab-seperated
                Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
                         17	        4639699	    A>T	    GLTPD2	n/a	    n/a	    34560
                         6	        25885022	T>C	    SLC17A4	n/a	    n/a	    55792
                         19	        49523688	T>A	    ZFP112	32218	39057	59962

    2. coverage: tab-separated-value file with the sequencing depth at the position of this variant in each sample
       DP (Total Depth) in VCF file (e.g. DP: 16)
       DP = 16

            coverage table: tab-seperated
                Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
                         17	        4639699	    A>T	    GLTPD2	n/a	    n/a	    34560
                         6	        25885022	T>C	    SLC17A4	n/a	    n/a	    55792
                         19	        49523688	T>A	    ZFP112	32218	39057	59962

    '''
    samp = [os.path.join(vcfdir, f) for f in os.listdir(vcfdir) if f.endswith('.vcf.gz')] #read all vcf files from folder
    df = pd.DataFrame() #final treeomics input files dataframe
    i = 1
    for s in samp:#extracting Allelic Depth (AD) and Total Depth (DP) from each vcf files and merging together
        print(f'Reading {s}...')
        d = pd.read_csv(s, sep="\t", header=0, comment = '##', engine='python', encoding = 'unicode_escape')
        d = d[d.FILTER == 'PASS'] #Extract variant with Filter: 'PASS'
        d = d[ ['#CHROM', 'POS', 'REF', 'ALT', 'INFO', 'FORMAT', d.columns[-1]] ]
        sname = d.columns[-1] #sample name
        d[f'{sname}AD'] = d[sname].str.split(':').str[1].str.split(',').str[1] #extract Allelic Depth (AD)
        d[f'{sname}DP'] = d[sname].str.split(':').str[3] #extract Total Depth (DP)
        d['Gene'] = d.INFO.str.split('[').str[1].str.split('|').str[0] #Extract Gene Symbol
        d['Change'] = d['REF'].astype(str) + '>' + d['ALT'].astype(str) #Merging Ref > Alt
            
        if(i == 1):#Only First VCF file Extracted information merge with empty df 
            df = d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']]
            i+=1
        else:#Merge other VCF file information
            df = df.merge(d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']], 
                          how='outer',
                          on=['#CHROM', 'POS', 'Change', 'Gene'])
           
    df.rename(columns= {'#CHROM': 'Chromosome', 'POS': 'Position'}, inplace=True)#change column names
    print(df[df.columns[~df.columns.str.contains('DP')]].head(4))#check output of merged VCF extracted information
    df[df.columns[~df.columns.str.contains('DP')]].to_csv(f"{odir}mut_read_table.txt", sep="\t", index=False)#remove columns with DP partial match and write into file
    df[df.columns[~df.columns.str.contains('AD')]].to_csv(f"{odir}coverage.txt", sep="\t", index=False)#remove col with AD partial match and write into file
    
    print(f"Successfully generated:\n\t1: {odir}mut_read_table.txt\n\t2: {odir}coverage.txt")

def treeomics_input_vcffilter(vcfdir='vcf', odir='./'):
    '''
    It generate input files from the VCF files from multiple samples located in one directory. 
    
    COMMAND:
    treeomics_input(vcfdir = './vcf/, odir='./)

    INPUT: Directory with all Funcotator Annotated GZ VCF Files

    OUTPUT:
    Generate Two Input files for TREEOMICS from Annotated VCF file

    1. <mut_reads.txt> mut-reads: number of reads reporting a variant (row) in each sample
       AD (Allele Depth) in VCF file (e.g. AD: 1,15; 1: ref and 15: alt)
       AD = 15

            mut-reads table: tab-seperated
                Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
                         17	        4639699	    A>T	    GLTPD2	n/a	    n/a	    34560
                         6	        25885022	T>C	    SLC17A4	n/a	    n/a	    55792
                         19	        49523688	T>A	    ZFP112	32218	39057	59962

    2. coverage: tab-separated-value file with the sequencing depth at the position of this variant in each sample
       DP (Total Depth) in VCF file (e.g. DP: 16). In the Gene column 'gene' is used suggests: no annotation
       DP = 16

            coverage table: tab-seperated
                Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
                         17	        4639699	    A>T	    gene	n/a	    n/a	    34560
                         6	        25885022	T>C	    gene	n/a	    n/a	    55792
                         19	        49523688	T>A	    gene	32218	39057	59962

    '''
    samp = [os.path.join(vcfdir, f) for f in os.listdir(vcfdir) if f.endswith('.vcf.gz')] #read all vcf files from folder
    df = pd.DataFrame() #final treeomics input files dataframe
    i = 1
    for s in samp:#extracting Allelic Depth (AD) and Total Depth (DP) from each vcf files and merging together
        print(f'Reading {s}...')
        d = pd.read_csv(s, sep="\t", header=0, comment = '##', engine='python', encoding = 'unicode_escape')
        d = d[d.FILTER == 'PASS'] #Extract variant with Filter: 'PASS'
        d = d[ ['#CHROM', 'POS', 'REF', 'ALT', 'INFO', 'FORMAT', d.columns[-1]] ]
        sname = d.columns[-1] #sample name
        d[f'{sname}AD'] = d[sname].str.split(':').str[1].str.split(',').str[1] #extract Allelic Depth (AD)
        d[f'{sname}DP'] = d[sname].str.split(':').str[3] #extract Total Depth (DP)
        d['Gene'] = 'gene' #No gene annotation so 'gene' is used
        #d['Gene'] = d.INFO.str.split('[').str[1].str.split('|').str[0] #Extract Gene Symbol
        d['Change'] = d['REF'].astype(str) + '>' + d['ALT'].astype(str) #Merging Ref > Alt
            
        if(i == 1):#Only First VCF file Extracted information merge with empty df 
            df = d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']]
            i+=1
        else:#Merge other VCF file information
            df = df.merge(d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']], 
                          how='outer',
                          on=['#CHROM', 'POS', 'Change', 'Gene'])
           
    df.rename(columns= {'#CHROM': 'Chromosome', 'POS': 'Position'}, inplace=True)#change column names
    print(df[df.columns[~df.columns.str.contains('DP')]].head(4))#check output of merged VCF extracted information
    df[df.columns[~df.columns.str.contains('DP')]].to_csv(f"{odir}mut_read_table.txt", sep="\t", index=False)#remove columns with DP partial match and write into file
    df[df.columns[~df.columns.str.contains('AD')]].to_csv(f"{odir}coverage.txt", sep="\t", index=False)#remove col with AD partial match and write into file
    
    print(f"Successfully generated:\n\t1: {odir}mut_read_table.txt\n\t2: {odir}coverage.txt")
    print(f'Before running treeomics_run:\n\tChange columns with sample name by one underscore, e.g. P1_0, P2_1 etc.')
