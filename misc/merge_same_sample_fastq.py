import os

def merge_same_sample(fastq_dir):
    '''
    INPUT: 'dir/path/fastq/files'
    OUTPUT: 'output/to/input/dir'
    RUN:
        merge_same_sample('dir/path/fastq/files')
    DESCRIPTION:
        if one sample has multiple fastq.gz files then it merge and provide file with name _merge_R1_ and delete the used fastq files
        e.g. A/a1_samp1_R1_001.fastq.gz
             A/a1_samp2_R1_001.fastq.gz
        output: A/a1_merge_R1_001.fastq.gz
    '''
    samples = [s.split('_')[0] for s in os.listdir(fastq_dir)]
    
    for s in set(samples):
        files_to_merge = [filename for filename in os.listdir(fastq_dir) if s in filename]
        
        # create empty lists to store the R1 and R2 files
        R1_files = [filename for filename in files_to_merge if 'R1' in filename]
        R2_files = [filename for filename in files_to_merge if 'R2' in filename]
        if len(R1_files) >1:
            os.system('cat ' + ' '.join([fastq_dir + '/' + x for x in R1_files])+' >' + fastq_dir + '/' + s + '_merge_R1_001.fastq.gz')
            print('cat ' + ' '.join([fastq_dir + '/' + x for x in R1_files])+' >' + fastq_dir + '/' + s + '_merge_R1_001.fastq.gz')
            os.system('rm ' + ' '.join([fastq_dir + '/' + x for x in R1_files]))
        if len(R2_files) >1:
            os.system('cat ' + ' '.join([fastq_dir + '/' + x for x in R2_files])+' >' + fastq_dir + '/' + s + '_merge_R2_001.fastq.gz')
            print('cat ' + ' '.join([fastq_dir + '/' + x for x in R2_files])+' >' + fastq_dir + '/' + s + '_merge_R2_001.fastq.gz')
            os.system('rm ' + ' '.join([fastq_dir + '/' + x for x in R2_files]))

        print(f'{s} is successfully merged...')
        

merge_same_sample('fastq')

