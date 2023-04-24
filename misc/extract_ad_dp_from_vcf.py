import os, pandas as pd, gzip

def extract_ad_pd(idir,odir,chr,start,end,fout):
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
        {odir}{fout}.txt #AD, AF, DP of SNPs in user defined genomic regions e.g. chr17, start, end
    RUN:
        extract_ad_pd(idir='filter/',odir='./',chr='chr17',start=7569720,end=7592868)
    '''
    files = [f for f in os.listdir(idir) if f.endswith('.filter.vcf.gz')]
    table = pd.DataFrame()
    for f in files:
        cols = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',f.split('.')[0]]
        df = pd.read_csv(f'{idir}{f}', sep='\t', comment='#', header=None, names=cols)
        df = df.assign(AD=df[f.split('.')[0]].str.split(':').str[1], 
                       AF=df[f.split('.')[0]].str.split(':').str[2], 
                       DP=df[f.split('.')[0]].str.split(':').str[3])
        df.rename(columns={'AD': 'AD_'+f.split('.')[0],
                           'AF': 'AF_'+f.split('.')[0],
                           'DP': 'DP_'+f.split('.')[0]}, inplace=True)
        df = df.filter(regex='CHROM|POS|AD|AF|DP')
        
        if table.empty: table = df
        else: table = pd.merge(table, df, on=['CHROM', 'POS'], how='outer')#; print(table);exit()
    table.to_csv(f'{odir}all_snps.txt',header=True,index=False)
    
    condition = (table['CHROM'] == chr[3:]) & (table['POS'] >= start) & (table['POS'] <= end)
    table[condition].to_csv(f'{odir}{fout}',sep="\t",header=True,index=False)
    print(f'files generated: \n\t{odir}{fout}\n\t{odir}all_snps.txt')

    
#extract_ad_pd(idir='filter/',odir='./',chr='chr17',start=7569720,end=7592868,fout='table.txt')