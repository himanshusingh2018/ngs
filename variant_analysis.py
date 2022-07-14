import os
from gatk_processing.vcfannotation import  anotfuncotator

def main():
    '''GATK VARIANT ANNOTATION'''
    #Variant Annotation using funcotator
    anotfuncotator(idir, odir, hgref, hgver, funcotator_data, ivcf, ovcf)
    print(f"{odir}{ovcf} is generated successfully...")
 

#data resources
idir='/Volumes/lilac_data_ziv/transciptome/paired_pnet/Project_12502_F/vcfpass/'
odir='./'
hgref='/Volumes/lilac_home_singhh5/hg38/grch38_gencode_bwa_index/GRCh38.primary_assembly.genome.fa'
hgver='hg38'
funcotator_data='/Volumes/lilac_data_ziv/transciptome/paired_pnet/funcotator_dataSources.v1.7.20200521s'
ivcf = 'M1.vcf' 
ovcf='1.vcf'

if __name__ == '__main__':
    main()