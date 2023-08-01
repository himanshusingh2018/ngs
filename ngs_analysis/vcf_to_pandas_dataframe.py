import pandas as pd, gzip

def vcf_to_pandas_dataframe(fvcf):
    '''
    Read VCF file and return into pandas dataframe
    INPUT: fvcf with path. e.g. vcf/P-0001993.vcf.gz
    '''
    if fvcf.endswith('.gz'):
        header = next(line for line in gzip.open(fvcf, 'rt', encoding='utf-8') if line.startswith('#CHROM')).strip().split('\t')
    else:
        header = next(line.strip() for line in open(fvcf, 'r', encoding='utf-8-sig') if line.startswith('#CHROM')).strip().split('\t')

    df = pd.read_csv(fvcf, sep='\t', comment='#', names=header, encoding='utf-8-sig')
    return(df)
#print(vcf_to_pandas_dataframe('vcf/P-0001993.vcf.gz'))

