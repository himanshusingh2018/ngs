import os

def bwa_fastq_aln(idir, odir, genome_index, reads, sample, nthreads):
    '''
    BWA MEM FastQ Alignment (preferred for variant )
        
        Variant Analysis Pipeline: https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/bwa-aligner.html
        Read Groups: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
            RG: Read Group; ID: Read Group Identifier; PU: Platform Unit; SM: Sample; PL: Platform/technoloty used to produce the read; LB: DNA preparation library

        BWA MEM Alignment Command:
            bwa mem -t {10} -R {RG} {genome_index} {r1.fastq (SE) | r1.fq r2.fq (PE)} | samtools view -bS --threads 10| samtools sort -o {/path/output/dir/r.sort.bam --threads 10}
                mem: 
                -t : number of threads
                -R : read group header line such as '@RG\tID:foo\tSM:bar'
            samtools:
                view: print all alignments if with no options
                -b  : Output BAM
                 S  : Ignored (input format is auto-detected)
                sort: Sort alignments by leftmost coordinates
                -o  : output e.g. sample.sort.bam
    Run:
        bwa_fastq_aln(genome_index = /path/to/genome_index/index.fa, reads = ['r1.fq', 'r2.fq], 
                        odir = '/output/path/', sample = 'sample_name', nthreads = 10)

    Arguments:
        genome_index: /path/to/bwa/genome/index_name
        reads       : list of reads. e.g. ['r1.fq'] or ['r1.fq', 'r2.fq']
        odir        : /output/dir/path/
        sample      : sample_name e.g. 'r'
        nthreads    :  10
    
    Output:
        /output/dir/path/sample.sort.bam
    '''
    os.makedirs(odir, exist_ok=True)#create output directory
    
    if(len(reads) == 1):
        RG = f"'@RG\\tID:{sample}_SE\\tPL:Illumina\\tPU:{sample}\\tLB:{sample}\\tSM:{sample}'"#ID: SE, bam header
        os.system(f"bwa mem -t {nthreads} -R {RG} {genome_index} {idir}{reads[0]} | samtools view -bS --threads {nthreads} | samtools sort -o {odir}{sample}.sort.bam --threads {nthreads}")
        print('Single-End alignment:\n\t{odir}{sample}.sort.bam')
    else:
        RG = f"'@RG\\tID:{sample}_PE\\tPL:Illumina\\tPU:{sample}\\tLB:{sample}\\tSM:{sample}'"#ID: PE, bam header
        os.system(f"bwa mem -t {nthreads} -R {RG} {genome_index} {idir}{reads[0]} {idir}{reads[1]} | samtools view -bS | samtools sort -o {odir}{sample}.sort.bam")
        print('Paired-End alignment:\n\t{odir}{sample}.sort.bam')
    
#bwa_fastq_aln(genome_index='/Volumes/lilac_home_singhh5/hg38/grch38_gencode_bwa_index/GRCh38.primary_assembly.genome.fa', reads=['fastq/EZ-123-B_IGO_08138_J_2_S101_R1_001.fastq.gz'], odir='./', sample='E', nthreads=4)

def bowtie2_fastq_aln(idir, odir, genome_index, reads, sample, nthreads):
    '''
    bowtie2: FastQ Alignemnt (prefer other than Variant Calling)
    https://rnnh.github.io/bioinfo-notebook/docs/bowtie2.html
    
        bowtie2 --no-unal --threads 4 -x genome_index_name -1 reads_1.fastq -2 reads_2.fastq | samtools view -bS --threads 4 | samtools sort -o {/path/outdir/}{sample}.sort.bam --threads 4
            --no-unal   : (optional) reads that do not align to the reference genome will not be written to sam output
            --threads   : number of threads
            -x          : genome index
            -1          : read1.fq
            -2          : read2.fq
    
    RUN: 
        bowtie2_fastq_aln(idir = '/path/bam/dir/', odir='/path/outdir/', genome_index = '/bowtie2/index/name', reads['r1.fq','r2.fq'], sample='sample')
    
    Arguments:
            idir            : /dir/path/input/bam
            odir            : /outdir/path/
            genome_index    : /path/bowtie2/index/GRCh38, bowtie2 genome index name 'GRCh38'
            reads           : list of read sequence, {['r1.fq'] or ['r1.fq', 'r2.fq']} 
            sameple         : name of sample
    
    Output:
        /outdir/path/sample.sort.bam
    '''
    os.makedirs(odir, exist_ok=True)#create output directory
    if(len(reads) == 1):
        RG = f"'@RG\\tID:{sample}_SE\\tPL:Illumina\\tPU:{sample}\\tLB:{sample}\\tSM:{sample}'"#ID: SE, bam header
        os.system("bowtie2 --no-unal --threads {nthreads} --rg {RG} -x {genome_index} {idir}{reads[0]} | samtools view -bS --threads {nthreads} | samtools sort -o {odir}{sample}.sort.bam --threads {nthreads}")
        print('Single-End alignment:\n\t{odir}{sample}.sort.bam')
    else:
        RG = f"'@RG\\tID:{sample}_SE\\tPL:Illumina\\tPU:{sample}\\tLB:{sample}\\tSM:{sample}'"#ID: SE, bam header
        os.system("bowtie2 --no-unal --threads {nthreads} --rg {RG} -x {genome_index} -1 {idir}{reads[0]} -2 {idir}{reads[1]} | samtools view -bS --threads {nthreads} | samtools sort -o {odir}{sample}.sort.bam --threads {nthreads}")
        print('Paired-End alignment:\n\t{odir}{sample}.sort.bam')

#bowtie2_fastq_aln(idir='./', odir='./', genome_index='GRCh38_noalt_as/GRCh38_noalt_as', reads=['fastq/EZ-123-B_IGO_08138_J_2_S101_R1_001.fastq.gz'], sample='E', nthreads=4)    

def star_rnaseq_fastq_2pass_aln(idir, odir, reads, genome_index, gtf, sample, nthreads):
    '''
    STAR 2 PASS FastQ Alignment
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/
    
    RUN: star_rnaseq_fastq_2pass_aln(idir='./', odir='./', reads=['f1.fq', 'f2.fq'], genome_index = '/path/to/genome_index', gtf = '/path/to/genome.gtf', reads=['r1.fq', 'r2.fq'], sample='E', nthreads=4)
                --runThreadN        : number of threads
                --genomeDir         : /path/to/genome_index
                --sjdbGTFfile       : /path/to/gtf_file
                --sjdbOverhang      : default: 100, 
                --readFilesIn       : specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. for Illumina 2x100b paired-end reads, the ideal value is 100-1=99 
                --readFilesCommand  : zcat. uncompressionCommand, sends uncompressed output to stdout
                --outFileNamePrefix : {odir}{sample}_SE_; output file prefix
                --outSAMtype        : BAM SortedByCoordinate Unsorted. output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command
                --twopassMode Basic : 1st pass mapping: it will automatically extract junctions, insert them into the genome index
                                      2nd pass mapping: re-map all reads

    '''
    os.makedirs(odir, exist_ok=True)#create temporary output directory 
    if(len(reads) == 1):
        os.system(f"STAR --runThreadN {nthreads} --genomeDir {genome_index} --sjdbGTFfile {gtf} --sjdbOverhang 100 --readFilesIn {reads[0]} --readFilesCommand zcat --outFileNamePrefix {odir}{sample}_SE_ --outSAMtype BAM SortedByCoordinate Unsorted --twopassMode Basic")
        os.system(f"samtools index {odir}{sample}_SE_Aligned.sortedByCoord.out.bam")
    else:
        os.system(f"STAR --runThreadN {nthreads} --genomeDir {genome_index} --sjdbGTFfile {gtf} --sjdbOverhang 100 --readFilesIn {reads[0]} {reads[1]} --readFilesCommand zcat --outFileNamePrefix {fout}{sample}_PE_ --outSAMtype BAM SortedByCoordinate Unsorted --twopassMode Basic")
        os.system(f"samtools index {odir}{sample}_PE_Aligned.sortedByCoord.out.bam")
