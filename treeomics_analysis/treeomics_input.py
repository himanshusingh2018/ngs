import os
import pandas as pd

def treeomics_input_by_vcf(idir, odir, featCoord, fvcf, filter=None):
    '''
    It generate input files from the VCF file
    
    COMMAND:
    treeomics_input_by_vcf(idir='./', odir='./', fvcf=fvcf, filter=None)
        
        idir   : /path/input/dir/
        odir   : /path/output/dir/
        fvcf   : 'file.vcf.gz'
        filter : None (not filtered); 'PASS': (filter the pass SNPs)

    INPUT: VCF File with all Samples, Mutect2 Variant Calling

    OUTPUT: Generate Two Input files for TREEOMICS from Annotated VCF file

    1. <mut_reads.txt> mut-reads: number of reads reporting a variant (row) in each sample
       AD (Allele Depth) in VCF file (e.g. AD: 1,15; 1: ref and 15: alt)
       AD = 15
                mut-reads table: tab-seperated
                    Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
                             17	        4639699	    A>T	    GLTPD2	n/a	    n/a	    34560
                             6	        25885022	T>C	    SLC17A4	n/a	    n/a	    55792

    2. coverage: tab-separated-value file with the sequencing depth at the position of this variant in each sample
       DP (Total Depth) in VCF file (e.g. DP: 16)
       DP = 16
                coverage table: tab-seperated
                    Columns: Chromosome Position    Change  Gene    Sample1 Sample2 Sample3
                             17	        4639699	    A>T	    GLTPD2	n/a	    n/a	    34560
                             6	        25885022	T>C	    SLC17A4	n/a	    n/a	    55792
    '''
    
        featCoord = pd.read_csv(featCoord, sep=",", header=0)
    d = pd.read_csv(f'{idir}{fvcf}', sep="\t", header=0, comment = '##', engine='python', encoding = 'unicode_escape')
    d = d if filter==None else d[d.FILTER=='PASS']#extract on the basis of filter None: extract all; other val: extract PASS Variants

    for sname in d.columns[9:d.shape[1]]:
        d[f'{sname}AD'] = d[sname].str.split(':').str[1].str.split(',').str[1] #extract Allelic Depth (AD)
        d[f'{sname}DP'] = d[sname].str.split(':').str[3] #extract Total Depth (DP)

    mapping_dict = dict(zip(featCoord['chr'], featCoord['symbol']))#Create a mapping dictionary from featCoord DataFrame

    def check_range(row):# Function to check if POS falls within start and end range
        chrom = row['#CHROM']
        pos = row['POS']
        matching_rows = featCoord[(featCoord['chr'] == chrom) & (featCoord['start'] <= pos) & (featCoord['end'] >= pos)]
        if len(matching_rows) > 0:
            return matching_rows.iloc[0]['symbol']
        return ''

    d['Gene'] = d.apply(check_range, axis=1)# Apply the check_range function row-wise to create the 'Gene' column
    d.rename(columns= {'#CHROM': 'Chromosome', 'POS': 'Position'}, inplace=True)#change column names
    d = d.explode('Gene')#convert list elements into rows
    d['Change'] = d['REF'].astype(str) + '>' + d['ALT'].astype(str) #Merging Ref > Alt
    d = d.assign(ALT=d.ALT.str.split(',')).explode('ALT').reset_index(drop=True)#all comma separated values in ALT in new rows
    d.rename(columns= {'#CHROM': 'Chromosome', 'POS': 'Position'}, inplace=True)#change column names
    #wiret coverage.txt and mut_read_table.txt
    d[['Chromosome','Position','Change','Gene'] + [x for x in d.columns if 'AD' in x]].to_csv(f"{odir}mut_read_table.txt", sep="\t", index=False)#remove columns with DP partial match and write into file
    d[['Chromosome','Position','Change','Gene'] + [x for x in d.columns if 'DP' in x]].to_csv(f"{odir}coverage.txt", sep="\t", index=False)#remove columns with DP partial match and write into file
    
    print(f"Successfully generated:\n\t1: {odir}mut_read_table.txt\n\t2: {odir}coverage.txt")
    
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
    samp = [os.path.join(vcfdir, f) for f in os.listdir(vcfdir) if f.endswith('.vcf')] #read all vcf files from folder
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
    samp = [os.path.join(vcfdir, f) for f in os.listdir(vcfdir) if f.endswith('.vcf')] #read all vcf files from folder
    df = pd.DataFrame() #final treeomics input files dataframe
    i = 1
    for s in samp:#extracting Allelic Depth (AD) and Total Depth (DP) from each vcf files and merging together
        print(f'Reading {s}...')
        d = pd.read_csv(s, sep="\t", header=0, comment = '##', engine='python', encoding = 'unicode_escape')
        d = d[ ['#CHROM', 'POS', 'REF', 'ALT', 'INFO', 'FORMAT', d.columns[-1]] ]
        sname = d.columns[-1] #sample name
        d[f'{sname}AD'] = d[sname].str.split(':').str[1].str.split(',').str[1] #extract Allelic Depth (AD)
        d[f'{sname}DP'] = d[sname].str.split(':').str[3] #extract Total Depth (DP)
        d = d.assign(ALT=d.ALT.str.split(',')).explode('ALT').reset_index(drop=True)#all comma separated values in ALT in new rows

        d['Gene'] = 'gene' #No gene annotation so 'gene' is used
        d['Change'] = d['REF'].astype(str) + '>' + d['ALT'].astype(str) #Merging Ref > Alt
        
        if(i == 1):#Only First VCF file Extracted information merge with empty df 
            df = d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']]; i+=1
        else:#Merge other VCF file information
            df = df.merge(d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']], 
                          how='outer', on=['#CHROM', 'POS', 'Change', 'Gene'])

    df.rename(columns= {'#CHROM': 'Chromosome', 'POS': 'Position'}, inplace=True)#change column names
    print(df[df.columns[~df.columns.str.contains('DP')]].head(4))#check output of merged VCF extracted information
    df[df.columns[~df.columns.str.contains('DP')]].to_csv(f"{odir}mut_read_table.txt", sep="\t", index=False)#remove columns with DP partial match and write into file
    df[df.columns[~df.columns.str.contains('AD')]].to_csv(f"{odir}coverage.txt", sep="\t", index=False)#remove col with AD partial match and write into file
    
    print(f"Successfully generated:\n\t1: {odir}mut_read_table.txt\n\t2: {odir}coverage.txt")
    print(f'Before running treeomics_run:\n\tChange columns with sample name by one underscore, e.g. P1_0, P2_1 etc.')

def treeomics_input_file_gtfannotation(idir, odir, fmut, fcov, feat_coord):
    '''
    Annotation of VCF file Using gencode grch38
    '''
    mut = pd.read_csv(f'{idir}{fmut}', sep="\t", header=0)
    cov = pd.read_csv(f'{idir}{fcov}', sep="\t", header=0)
    feat_coord = pd.read_csv(feat_coord, sep=",", header=0)

    def gen_search(x): #extract gene symbol for the SNP Position of VCF file
        chrom, pos = x['Chromosome'], x['Position'] 
        return feat_coord[(feat_coord.chr == chrom) & (pos >= feat_coord.start) & (pos <= feat_coord.end)].symbol.values
    
    mut['Gene'] = pd.DataFrame({'symbol': mut[['Chromosome', 'Position']].head(1000).apply(gen_search, axis=1)})#add new column generated by gen_search
    mut = mut.explode('Gene')#convert list elements into rows
    
    #mut = mut.astype(int)
    mut.to_csv(f'{odir}{fmut[:-4]}.anot.txt', header=True, sep="\t", index=False)#write to the file
    
    cov['Gene'] = pd.DataFrame({'symbol': cov[['Chromosome', 'Position']].head(1000).apply(gen_search, axis=1)})#add new column generated by gen_search
    cov = cov.explode('Gene')#convert list elements into rows
    cov.to_csv(f'{odir}{fcov[:-4]}.anot.txt', header=True, sep="\t", index=False)#write to the file
    print(f'Successfully generated:\n\t{odir}{fmut[:-4]}.anot.txt\n\t{odir}{fcov[:-4]}.anot.txt')

def treeomics_input_vcfanotgtf(vcfdir='vcf', odir='./'):
    '''
    It generate input files from the VCF files from multiple samples located in one directory. 
    
    COMMAND:
    treeomics_input(vcfdir = './vcf/, odir='./')

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
    samp = [os.path.join(vcfdir, f) for f in os.listdir(vcfdir) if f.endswith('.anot.vcf')] #read all vcf files from folder
    df = pd.DataFrame() #final treeomics input files dataframe
    i = 1
    for s in samp:#extracting Allelic Depth (AD) and Total Depth (DP) from each vcf files and merging together
        print(f'Reading {s}...')
        d = pd.read_csv(s, sep="\t", header=0, comment = '##', engine='python', encoding = 'unicode_escape')
        #print(d.head());exit()
        #d = d[d.FILTER == 'PASS'] #Extract variant with Filter: 'PASS'
        d = d[ ['#CHROM', 'POS', 'REF', 'ALT', 'INFO', 'FORMAT', d.columns[-2], 'gene'] ]
        sname = d.columns[-2] #sample name
        d[f'{sname}AD'] = d[sname].str.split(':').str[1].str.split(',').str[1] #extract Allelic Depth (AD)
        d[f'{sname}DP'] = d[sname].str.split(':').str[3] #extract Total Depth (DP)
        #d['Gene'] = d.INFO.str.split('[').str[1].str.split('|').str[0] #Extract Gene Symbol
        d['Gene'] = d.gene #Extract Gene Symbol
        d['Change'] = d['REF'].astype(str) + '>' + d['ALT'].astype(str) #Merging Ref > Alt
            
        if(i == 1):#Only First VCF file Extracted information merge with empty df 
            df = d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']]
            i+=1
        else:#Merge other VCF file information
            df = df.merge(d[['#CHROM', 'POS', 'Change', 'Gene', f'{sname}AD', f'{sname}DP']], 
                          how='outer',
                          on=['#CHROM', 'POS', 'Change', 'Gene'])
           
    df = df.sort_values(by = ['#CHROM', 'POS'], ascending = [True, True])#sort by chromosome and position
    df.rename(columns= {'#CHROM': 'Chromosome', 'POS': 'Position'}, inplace=True)#change column names
    print(df[df.columns[~df.columns.str.contains('DP')]].head(4))#check output of merged VCF extracted information
    df[df.columns[~df.columns.str.contains('DP')]].to_csv(f"{odir}mut_read_table.txt", sep="\t", index=False)#remove columns with DP partial match and write into file
    df[df.columns[~df.columns.str.contains('AD')]].to_csv(f"{odir}coverage.txt", sep="\t", index=False)#remove col with AD partial match and write into file
    
    print(f"Successfully generated:\n\t1: {odir}mut_read_table.txt\n\t2: {odir}coverage.txt")
