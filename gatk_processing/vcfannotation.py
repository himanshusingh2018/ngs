import os
def download_funcotator_data(dtype='somatic'):
    '''
    Download latest version of funcotator_dataSources.v1.7.20200521s
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


def anotfuncotator(idir, odir, hgref, hgver, funcotator_data, ivcf, ovcf):
    '''
    VCF Annotation using Funcotator
    
    Arguments:
        idir: input directory; odir: output directory
        ivcf: input vcf file; ovcf: output vcf file
        reference: GRCh38 genome fasta file
        ref-version: hg38, hg19
        data-sources-path: funcotator_dataSources.v1.7.20200521s
    
    Run:
        gatk Funcotator --variant idir/in.vcf --reference GRCh38.fa --ref-version hg38 --data-sources-path funcotator_dataSources.v1.7.20200521s --output odir/out.vcf --output-file-format VCF
    '''
    os.system(f"gatk Funcotator --variant {idir}{ivcf} --reference {hgref} --ref-version {hgver} --data-sources-path {funcotator_data} --output {odir}{ovcf} --output-file-format VCF")
    print(f"{odir}{ovcf} generated successfully...")


def vcfTotable(idir, odir, ivcf, otable):
    '''
    VCF to Table
    Columns: CHROM  POS  REF  ALT  FUNCOTATION  AD  DP
        CHROM - Chromosome Number
        POS   - Position
        REF   - Reference Allele
        ALT   - Alternative Allele
        AD    - Allele Depth (e.g. AD: 1,15; 1: ref and 15: alt)
        DP    - Total Depth (e.g. DP: 16)

    Arguments:
        idir: input directory; odir: output directory
        ivcf: input vcf file; ovcf: output vcf file
    
    Run:
        gatk VariantToTable -V idir/in.vcf -F CHROM -F POS -F REF -F ALT -F FUNCOTATION -GF AD -GF DP -O odir/out.table
    '''
    os.system(f"gatk VariantsToTable -V {idir}{ivcf} -F CHROM -F POS -F REF -F ALT -F FUNCOTATION -GF AD -GF DP -O {odir}{otable}")
    print(f"{odir}{otable} generated successfully...")

anotfuncotator(idir='/Volumes/lilac_data_ziv/transciptome/paired_pnet/Project_12502_F/vcfpass/', 
                odir='./', 
                hgref='/Volumes/lilac_home_singhh5/hg38/grch38_gencode_bwa_index/GRCh38.primary_assembly.genome.fa', 
                hgver='hg38', 
                funcotator_data='/Volumes/lilac_data_ziv/transciptome/paired_pnet/funcotator_dataSources.v1.7.20200521s', 
                ivcf = 'M1.vcf', 
                ovcf='1.vcf')
