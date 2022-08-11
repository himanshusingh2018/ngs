import pandas as pd
def gtf_snp_annotation(genomefeature_coordinates, fvcf):
    '''
    Annotating SNPs
    '''
    genomefeat = pd.read_csv(genomefeature_coordinates, sep=",", header=0, names=['chr','start','end','strand','gene'])
    
    print(genomefeat)
    vcf = pd.read_csv(fvcf, sep="\t", header=0, comment = '##', engine='python', encoding = 'unicode_escape')
    print(vcf[['POS','#CHROM']])
    print(vcf.columns)

fvcf= '/Volumes/lilac_data_ziv/transciptome/paired_pnet/Project_12502_F/test_out/vcf/P-0032461-N01-WES_IGO_12502_F_15_S14_L001.vcf.gz'
genomefeature_coordinates='/Volumes/lilac_home_singhh5/software/out.3.txt'
gtf_snp_annotation(genomefeature_coordinates=genomefeature_coordinates, fvcf=fvcf)
