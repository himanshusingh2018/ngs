import pandas as pd
def gene_length(odir,gtf_url, feature):
    '''
    Extract Gene Length from UCSC GTF file
    Input: 
        gtf_url=hg19.refGene.gtf.gz
                http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.refGene.gtf.gz
        feature: transcript (gene length = intron+exon according to the transcript)
                    exon (all the exon of corresponding genes will be added)
        
        odir: /output/dir/
    Output:
        odir/{feature}.genelength.txt
        e.g. odir/exon.genelength.txt or odir/transcript.genelength.txt
    Run:
    import gene_length
    gene_length(idir = /path/hg19/gtf/
                odir = /path/to/output/
                gtf_url = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.refGene.gtf.gz',
                feature = 'exon' or 'transcript')
    '''
    gtf = pd.read_csv(gtf_url, sep='\t', compression='gzip', comment='#', header=None, names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    gtf['gene_symbol'] = gtf['attribute'].apply(lambda x: x.split(';')[0].split('"')[1])
    gtf['gene_length'] = (gtf.start - gtf.end).abs()+1
    
    if feature == 'transcript':
        gtf.drop_duplicates(subset='gene_symbol', keep='first', inplace=True)
        gtf.to_csv(f'{odir}{feature}.genelength.txt')
        print(f'File generated: {odir}{feature}.genelength.txt')
    else:
        gtf = gtf[gtf.feature=='exon'].groupby('gene_symbol')['gene_length'].agg(['sum', 'count'])
        gtf.rename(columns={'sum': 'total_exon_length', 'count': 'no_of_exons'}, inplace=True)
        gtf.to_csv(f'{odir}exon.genelength.txt')
        print(f'File generated: {odir}exon.genelength.txt')

    return(gtf)

#a = gene_length(odir='~/Desktop/',gtf_url='http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/genes/hg19.refGene.gtf.gz', feature='exon')