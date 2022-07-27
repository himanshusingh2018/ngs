import os
import pandas as pd

from gatk_processing.vcfannotation import  anotfuncotator
from treeomics_analysis.treeomics_input import treeomics_input
def main():
    '''GATK VARIANT ANNOTATION'''
    #Variant Annotation using funcotator
    #anotfuncotator(idir, odir, hgref, hgver, funcotator_data, ivcf, ovcf)
    #generate treeomics input files: 1) coverage.txt; 2) mut_read_table.txt
    treeomics_input(vcfdir = idir, odir=odir)
    print(f"{odir}coverage.txt is generated successfully...")
 

#data resources
idir='/Volumes/lilac_data_ziv/transciptome/paired_pnet/Project_12502_F/anotvcf/'
odir='./'
hgref='/Volumes/lilac_home_singhh5/hg38/grch38_gencode_bwa_index/GRCh38.primary_assembly.genome.fa'
hgver='hg38'
funcotator_data='/Volumes/lilac_data_ziv/transciptome/paired_pnet/funcotator_dataSources.v1.7.20200521s'


if __name__ == '__main__':
    main()