import pandas as pd
def gene_length(gtf_url):
    #gtf_url = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.refGene.gtf.gz'
    #gtf_url = '/Users/singhh5/Downloads/hg19.refGene.gtf.gz'
    gtf = pd.read_csv(gtf_url, sep='\t', compression='gzip', comment='#', header=None, names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    gtf['gene_symbol'] = gtf['attribute'].apply(lambda x: x.split(';')[0].split('"')[1])
    gtf['gene_length'] = (gtf.start - gtf.end).abs()+1
    gtf.drop_duplicates(subset='gene_symbol', keep='first', inplace=True)
    return(gtf)

a = gene_length(gtf_url='http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.refGene.gtf.gz')
print(a)