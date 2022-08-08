import os, sys, glob
import pandas as pd
import functools
from functools import reduce

def main():
    from ngs.ngs_analysis import picard_process as pc

    #STEP1: Picard Clean BAM files
    pc.picard_bam_clean(idir = bam_aln, odir= bam_mkdup, sample=f'{sname}.sort.bam')  

    #STEP2: Picard Mark Duplicates in Clean BAM files
    pc.picard_markduplicate(idir=bam_aln, odir=bam_mkdup, sample=f'{sname}.clean.bam')

    #STEP3: Copy BAM Files in Picard Duplicate Marked BAM Files
    print(f'cp {bam_mkdup}{sname}.mkdup.bam* {bam_aln}')

    #STEP4: Delete BAM if bam.mkdup is copied
    if os.path.exists(f'{bam_aln}{sname}.mkdup.bam'): os.remove(f'{bam_aln}{sname}.sort.bam*')



#data resources
bam_aln = '/data/ziv/transciptome/paired_pnet/Project_12502_F/bwa_aln/'
bam_mkdup = '/scratch/him/Project_12502_F/mkdup/'

index = sys.argv[1]
sname = [sample.split('.sort.bam')[0] for sample in os.listdir(bam_aln) if sample.endswith('.sort.bam')]

if __name__ == '__main__':
    if index == "TEST":
        print(f"set job index to [1-len(samples)]")
    else:
        main()

print('All processes are completed successfully...')