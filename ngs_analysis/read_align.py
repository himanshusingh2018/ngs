import os

def bwa_align(idir, fin, genome_size_in_gb, nthreads):
  '''
  BWA Genome Indexing (preferred for variant analysis)
  http://bio-bwa.sourceforge.net/bwa.shtml

  Run:
    bwa index -a is|bwtsw {idir=/path/to/}{fin = genome.fa} -p -t {nthreads = 5}
  
  Parameters:
    index: Index genome sequence in fasta format
    -a   : Alogorithm for constructing BWT index. Two algorithms bwtsw and is
    is   : Genome size < 2GB
    bwtsw: Genome size >= 2GB. 
    -p   : Prefix of the output database, same as genome file name. e.g. GRCh38.primary_assembly.genome.fa
    -t   : number of threads
    
  Arguments:
    idir : Genome fasta file directory path. e.g. /path/to/genome/
    fin  : Genome fasta file name. e.g. genome.fa
    genome_size_in_gb: 3 (Human genome size: 3.2 GB)
    nthreads: number of threads
  
  Output: Five output file will generated at the same location of genome fasta file with same name
          1. GRCh38.primary_assembly.genome.fa.amb
          2. GRCh38.primary_assembly.genome.fa.ann
          3. GRCh38.primary_assembly.genome.fa.bwt
          4. GRCh38.primary_assembly.genome.fa.pac
          5. GRCh38.primary_assembly.genome.fa.sa

  '''
  print(f'bowtie2 -t -x bowtie2/NC_012967.1 -1 SRR030257_1.fastq -2 SRR030257_2.fastq -S bowtie2/SRR030257.sam')


def bwa_fastq_aln(genome_index, read1, read2, fout, sample, nthreads):
    """ 
    STEP 2: BWA MEM FastaQ Alignment
    """    
    os.makedirs(fout, exist_ok=True)
    #fastq alignment
    '''
    bwa mem -t 24 -R "@RG\tID:P-0000879-N01-WES_IGO_12502_F_25_S21_L001" /home/singhh5/hg38/grch38_gencode_bwa_index/GRCh38.primary_assembly.genome.fa fastq/P-0000879-N01-WES_IGO_12502_F_25_S21_L001_R1_001.fastq.gz fastq/P-0000879-N01-WES_IGO_12502_F_25_S21_L001_R2_001.fastq.gz | samtools view -bS | samtools sort -o bwa_aln/P-0000879-N01-WES_IGO_12502_F_25_S21_L001.sort.bam
    '''
    #RG = "@RG\tID:{sample}_PE\tPL:Illumina\tPU:{sample}\tLB:{sample}\tSM:{sample}"#bam header
    RG = f"'@RG\\tID:{sample}_PE\\tPL:Illumina\\tPU:{sample}\\tLB:{sample}\\tSM:{sample}'"#bam header
    aln = f"bwa mem -t {nthreads} -R {RG} {genome_index} {read1} {read2} | samtools view -bS | samtools sort -o {fout}{sample}.sort.bam"
    print(aln)
    os.system(aln)
    exit()
    #indexing bam file
    os.system(f"samtools index {fout}{sample}.sort.bam")
    
    print(f'{sample} alignment and indexing is completed successfully...')

def picard_clean_n_markduplicate(sample, bwa_aln):
    '''
    PICARD CLEANING SORT BAM FILE
    setting MAPQ to 0 for unmapped reads
    '''
    os.system(f"picard CleanSam --INPUT {bwa_aln}{sample}.sort.bam --OUTPUT {bwa_aln}{sample}.clean.bam")
    print(f"{bwa_aln}{sample}.clean.bam is generated successfully...")
    os.system(f"picard MarkDuplicates --INPUT {bwa_aln}{sample}.clean.bam --OUTPUT {bwa_aln}{sample}.mkdup.bam --METRICS_FILE {bwa_aln}{sample}.dup.metrics.txt --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true")
    print(f"Successfully generated files: \n\t1: {bwa_aln}{sample}.sort.mkdup.bam \n\t2:{fout}{sample}.dup.metrics.txt")

def gatk_variant_calling_sample(sample, fout, nthreads):
    '''
    GATK Mutect2 Variant Calling
    '''
    os.system(f"gatk Mutect2 --native-pair-hmm-threads {nthreads} -R {genome_index} -I {fout}{sample}_R1_001.sort.bam -O {fout}{sample}.vcf.gz")
    print('{sample} is generated successfully:\n\t{fout}{sample}.vcf.gz')