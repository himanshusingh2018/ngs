import pandas as pd

vcf = pd.read_csv('~/Downloads/P-0032461-N01-WES_IGO_12502_F_15_S14_L001.vcf.gz', sep="\t", comment = '##', engine='python', encoding = 'unicode_escape')
genomefeat = pd.read_csv('~/Downloads/out.3.txt', sep=",", header=0)

print(vcf.head())
print(genomefeat.head())

def gen_search(x):
    chrom = x['#CHROM']
    pos = x['POS']
    return genomefeat[(genomefeat.chr == chrom) & (pos >= genomefeat.start) & (pos <= genomefeat.end)].symbol.values


vcf['gene'] = pd.DataFrame({'symbol': vcf[['#CHROM', 'POS']].head(1000).apply(gen_search, axis=1)})
vcf['gene'] = vcf['gene'].str.strip('[]').str.split(',\s*')
vcf['liststring'] = [','.join(map(str, l)) for l in vcf['gene']]

print(vcf)
vcf.to_csv('vcf.txt',header=True,index=False)











'''
def gen_search(x):
    chrom = x['#CHROM']
    pos = x['POS']
    return genomefeat[(genomefeat.chr == chrom) & (pos > genomefeat.start) & (pos < genomefeat.end)].symbol.values

df = pd.DataFrame({'symbol': vcf[['#CHROM', 'POS']].head(1000).apply(gen_search, axis=1)})
print(df.head())


# vcf['symbol'] = np.where(((vcf.POS >= genomefeat.start) | (vcf.POS <= genomefeat.end)) & (vcf['#CHROM'] == genomefeat.chr)

for f1 in vcf:
    for f2 in genomefeat:
        # if( f1["#CHROME"] == f2["chr"] and f1["POS"] >= f2["start"] and f1["POS"] <= f2["POS"] ):
        #     f1["GEN"]= f2["symbol"]
        pass

print(vcf.head())
print(genomefeat.head())

# vcf['symbol'] = np.where(((vcf.POS >= genomefeat.start) | (vcf.POS <= genomefeat.end)) & (vcf['#CHROM'] == genomefeat.chr)
'''