import os, sys, glob
import functools
def extract_data_between_two_pattern(idir, odir, fname, START, END, fout):
    '''
    Extract Data between two markers START and END

    RUN: 
    extract_data_between_two_pattern(idir='./',odir='./',
                                     fname = 'grch38/GRCh38.genome.fa',
                                     START='>chr21', END = '>chr22', 
                                     fout = 'grch38/GRCh38.chr21.fa')
    idir : /input/dir/path/
    odir : /output/dir/path/
    fname: input file name
    START: start pattern
    END  : end pattern
    fout : output file name
    '''
    os.system(f"sed -n '/{START}/,/{END}/p' {idir}{fname} | sed '$d' >{odir}{fout}")
    print(f'{odir}{fout} is generated successfully...')
    