import os, sys, glob
from functools import reduce

def htseq_exon_count(idir, odir, gtf, mode, sample):
    """
    STEP 4: HTSEQ Exon Count
    
    Mode to handle reads overlapping more than one feature
    choices: union, 
             intersection-strict, 
             intersection-nonempty
    """
    os.system(f'htseq-count -m {mode} -t exon -i gene_id -f bam {idir}{sample}.bam {gtf} >{odir}{sample}.union.htseq.count')
    print(f'HTSeq Exon Count of {odir}{sample}.union.htseq.count is successfully generated...')

def merge_all_htseqcount(idir, odir,fext):
    """
    STEP: Merge All HTSeq Count
    
    RUN: merge_all_htseqcount(idir='./', odir='./', fext='extension_of_file')

    Argument:
            idir: /dir/path/to/input/files/
            odir: /dir/path/to/output/files/
            fext: extension of input files
    Output:
        /output/dir/allSample.htseq.count.txt  
    """
    flist = glob.glob(os.path.join(idir, f"*.{fext}")) #fpath with file name
    df = reduce(lambda l, r: l.merge(r, on='ENSG'), 
                        [pd.read_csv(f,sep="\t",
                         header=None, 
                         names=['ENSG',f.split('/')[-1][:-18] ]) for f in flist ])
    df = df[df.ENSG.str.contains('ENSG')]
    df.to_csv(f'{odir}allSample.htseq.count.txt',sep="\t",header=True,index=False)