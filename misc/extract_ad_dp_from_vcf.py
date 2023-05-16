import os, gzip
import pandas as pd

def extract_ad_pd(idir, odir, chr, start, end, fout):
    '''
    Extract Allelic Depth and Depth from filtered VCFs in DIR
    INPUT:
        idir: /path/filter/vcf/dir/ #gatk4 filtered vcf file
        odir: /path/of/output/dir/
        chr: chrNo e.g. chr17
        start: start of genomic coordinate
        end: end of genomic coordinate
        fout: name of output file
    OUTPUT:
        {odir}all_snp.txt #AD, AF, DP of all SNP position in samples vcf in {idir}
        {odir}{fout}.txt #AD, AF, DP of SNPs in user-defined genomic regions e.g. chr17, start, end
    RUN:
        extract_ad_pd(idir='filter/', odir='./', chr='chr17', start=7569720, end=7592868, fout='table.txt')
    '''
    files = [f for f in os.listdir(idir) if f.endswith('.filter.vcf.gz')]
    table = pd.DataFrame()
    for f in files:
        cols = line = next((line.rstrip('\n').split('\t') for line in gzip.open(f'filter/{f}', 'rt', encoding='latin-1') if line.strip().startswith('#CHROM')), None)
        df = pd.read_csv(f'{idir}{f}', sep='\t', comment='#', header=None, names=cols)
        for i in range(9, len(cols)):#iterate over each sample in vcf file
            df = df.assign(
                AD=df[cols[i].split('.')[0]].str.split(':').str[1],
                AF=df[cols[i].split('.')[0]].str.split(':').str[2],
                DP=df[cols[i].split('.')[0]].str.split(':').str[3]
            )
            df.rename(
                columns={
                    'AD': 'AD_'+cols[i],
                    'AF': 'AF_'+cols[i],
                    'DP': 'DP_'+cols[i]
                },
                inplace=True
            )
        df = df.filter(regex='CHROM|POS|AD|AF|DP')
        
        df['POS'] = df['POS'].astype(int)  # Convert 'POS' column to int64 data type
        
        if table.empty:
            table = df
        else:
            table = pd.merge(table, df, on=['#CHROM', 'POS'], how='outer')

    table.to_csv(f'{odir}all_snps.txt', header=True, index=False)

    condition = (table['#CHROM'] == chr[3:]) & (table['POS'] >= start) & (table['POS'] <= end)
    table[condition].to_csv(f'{odir}{fout}', sep="\t", header=True, index=False)
    print(f'Files generated:\n\t{odir}{fout}\n\t{odir}all_snps.txt')


#extract_ad_pd(idir='filter/', odir='./', chr='chr17', start=7569720, end=7592868, fout='table.txt')
