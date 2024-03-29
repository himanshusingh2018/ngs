import os
import pandas as pd
def gencode_feature_coordinates(gtf, idir='./', odir='./', fout='out'):
    '''
    Extract FEATURE Specific Coordinates from Gencode GTF file (Region: PRI)
    
    RUN: gencode_feature_coordinates(idir='./', odir='./', fout = 'out.txt', gtf='gencode.gtf', feature) 

    Input:
            idir   : /path/to/gtf/
            odir   : /path/to/output/dir/
            fout   : out
            gtf    :  GTF file downloaded from Gencode (https://www.gencodegenes.org/human/)
                      GTF Format: https://www.gencodegenes.org/pages/data_format.html
            
            feature:  {'transcript' | 'exon' | 'gene' | 'start_codon' | 'stop_codon' | 'UTR' | 'CDS' | 'Selenocysteine'}

    
    Output: Returns Feature Coordinates and Gene Name Annotation
            Output file: 'out.feature.txt'
            chr,start,end,strand,symbol
            chr1,11869,12227,+,DDX11L1
            chr1,12613,12721,+,DDX11L1
    '''
    feature =  {1:'transcript', 2:'exon', 3:'gene', 4:'start_codon', 5:'stop_codon', 6:'UTR', 7:'CDS', 8:'Selenocysteine'}
    for k, v in feature.items(): print(f'"{k}": {v.upper()}')#print all the options 
    option = int(input('Type numeric option: '))
    if option in feature.keys():
        print('\nSelected feature:', feature[option])
        f = pd.read_csv(f'{idir}{gtf}', sep="\t", header=None, comment='#', usecols=[0,2,3,4,6,8], names=['chr','feature','start','end', 'strand', 'info'])#read specific columns from gtf file
        feat = f[f.feature == feature[option]]#[['chr','genomic_start','genomic_end', 'strand', 'info']]
        feat['symbol'] = feat['info'].str.split('gene_name "').str[1].str.split('";').str[0]
        #feat['info'] = feat['info'].str.split(';', n=3, expand=False).str[2].str.slice(start=12, stop=-1, step=1)#extract only gene symbol from the column value
        feat[['chr','start','end', 'strand', 'symbol']].to_csv(f'{odir}{fout}.{feature[option]}.txt', index=False)
    else:
        print('Please enter correct Option'); exit()
    print(f'Successfully generated:\n\t{odir}{fout}.{feature[option]}.txt')

def gtf_vcfsnp_annotation(idir, odir, genomefeature_coordinates, fvcf, fout):
    '''
    Annotating SNPs By Extracted misc.gtf_snp_annotation Module

    RUN:
        gtf_vcfsnp_annotation(idir = './', odir='./', 
                              genomefeature_coordinates='/path/misc.gtf_process.gencode_feature_coordinates/output/out.txt, 
                              fvcf = 'snp.vcf', fout='output.txt')

        idir: /path/input/vcf/
        odir: /path/outfile
        genomefeature_coordinates: /path/to/misc.gtf_process.gencode_feature_coordinates/output/out.txt
        fvcf: 'input file, mutect2.snp.vcf'
        fout: 'output file name, annotated gene symbol vcf file

            vcf.anot.txt: (Comma separated mutect2 vcf annotated file; One column gene (gene symbol) is added)
                #CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,P-0032461-N01-WES_IGO_12502_F_15_S14_L001,gene
                chr1,13613,.,T,A,.,.,"AS_SB_TABLE=8,31|3,0;DP=42;ECNT=1;MBQ=37,37;MFRL=276,305;MMQ=23,22;MPOS=34;POPAF=1.70;TLOD=4.19",GT:AD:AF:DP:F1R2:F2R1:FAD:SB,"0/1:39,3:0.090:42:18,1:19,2:39,3:8,31,3,0",DDX11L1
                chr1,13813,.,T,G,.,.,"AS_SB_TABLE=0,8|0,4;DP=12;ECNT=2;MBQ=37,37;MFRL=312,295;MMQ=22,22;MPOS=38;POPAF=1.66;TLOD=14.19",GT:AD:AF:DP:F1R2:F2R1:FAD:PGT:PID:PS:SB,"0|1:8,4:0.357:12:4,2:4,2:8,4:0|1:13813_T_G:13813:0,8,0,4",DDX11L1
    '''
    vcf = pd.read_csv(f'{idir}{fvcf}', sep="\t", comment = '##', engine='python', encoding = 'unicode_escape')
    genomefeat = pd.read_csv(genomefeature_coordinates, sep=",", header=0)

    def gen_search(x): #extract gene symbol for the SNP Position of VCF file
        chrom, pos = x['#CHROM'], x['POS'] 
        return genomefeat[(genomefeat.chr == chrom) & (pos >= genomefeat.start) & (pos <= genomefeat.end)].symbol.values

    vcf['gene'] = pd.DataFrame({'symbol': vcf[['#CHROM', 'POS']].head(1000).apply(gen_search, axis=1)})#add new column generated by gen_search
    vcf = vcf.explode('gene')#convert list elements into rows
    vcf.to_csv(f'{odir}{fout}', header=True, sep="\t", index=False)#write to the file
    print(f'Successfully generated:\n\t{odir}{fout}')