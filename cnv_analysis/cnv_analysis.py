import os
import pandas as pd


def accessible_region_gencode_ref_genome(idir, odir, fin, fout):
    '''
    First step is to create a bed file that contains the locations of the accessible regions for a given reference genome.
    Gencode is preferred for Genome GTF and FASTA

    Input:
        gencode genome fasta file: GRCh37.fa
        
    Run:
        cnvkit.py access /path/to/gencode.hg19.ref.fa -o /path/to/gencode.hg19.ref.accessible.bed

    Output: gencode.hg19.ref.accessible.bed
    '''

    os.system(f"cnvkit.py access {idir}{fin} -o {odir}{fout}.accessible.bed")
    print(f'Successfully generated:\n\t{odir}{fout}.accessible.bed')

def refflet_gtf(idir, odir, ingtf, fout):
    '''
    UCSC refFlat.txt download from google cloud
    https://cloud.google.com/life-sciences/docs/resources/public-datasets/ucsc?_ga=2.45345068.-980710568.1698679471
    '''
    os.system(f'wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz')
    #pass




def cnv_tumor_only():
    '''
    CNV Analysis
    Link: https://sequencing.roche.com/content/dam/rochesequence/worldwide/resources/SPR-White-Paper-KAPA-TE-data-evaluation-for-somatic-mutations.pdf
    
    Argument:
        dtype: somatic or germline

    Run:
        download_funcotator_data(dtype='somatic')
    '''
    if(dtype == 'somatic'):
        os.system(f"gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download")
        print("Funcotator somatic data is downloaded successfully...")
    else:
        os.system(f"gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download")
        print("Funcotator germline data is downloaded successfully...")
