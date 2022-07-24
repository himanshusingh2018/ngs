import os

def accessible_region_ensembl_ref_genome(idir, odir, fin, fout):
    '''
    First step is to create a bed file that contains the locations of the accessible regions for a given reference genome.
    
    Run:
        cnvkit.py access /path/to/ensembl.hg19.ref.fa -o /path/to/ensembl.hg19.ref.accessible.bed
    
    Output:

    '''
    os.system(f"cnvkit.py access {idir}{fin}.fa -o {odir}{fout}.accessible.bed")


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