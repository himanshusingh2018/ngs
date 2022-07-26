import os

def bwa_gi(idir, fin, genome_size_in_gb, nthreads):
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
  if(genome_size_in_gb < 2):
    os.system(f'bwa index -a is {idir}{fin} -p -t {nthreads}')
  else:
    os.system(f'bwa index -a bwtsw {idir}{fin} -p -t {nthreads}')
  print(f'BWA indexes generated successfully:\n\t{idir}{fin}')

def bowtie2_gi(idir, fin, name_index, nthreads):
  '''
  bowtie2 Genome Indexing
  https://rnnh.github.io/bioinfo-notebook/docs/bowtie2.html

  bowtie2 has no upper limit on read length, can make gapped alignments, more flexible for paired-end alignment, faster and more memory efficient
  bowtie is advantageous over bowtie2 for relatively short sequencing reads (50bp or less)
  
  RUN:
    bowtie2-build {idir=/path/to/}{fin=genome.fa} {idir=/path/to/}{name_index=GRCh38} --threads {nthreads = 10}
  
  Arguments:
    idir      : Genome file directory path
    fin       : Genome fasta file
    name_index: name prefix of the index files
    nthreads   : number of threads

  Output: Six output file will generated at the same location of genome fasta file with same name
          1. index_name.1.bt2
          2. index_name.2.bt2
          3. index_name.3.bt2
          4. index_name.4.bt2
          5. index_name.rev.1.bt2
          6. index_name.rev.2.bt2

  '''
  print(f'bowtie2-build {idir}{fin} {idir}{name_index} --threads {nthreads}')

def star_gi(idir, fin, nthreads):
  '''
  STAR Genome Index
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/

  RUN:
    STAR --runThreadN {nthreads=8} --runMode genomeGenerate --genomeDir {fin='./genome/'} --genomeFastaFiles {fin}{genome.fa}
  
  Arguments:
    idir      : Genome file directory path
    fin       : Genome fasta file
    name_index: name prefix of the index files
    nthreads  : number of threads

  Output: Eleven output file will be generated at the same location of genome fasta file with same name, ##ls -sh1
           1.  4.0K chrLength.txt
           2.  4.0K chrNameLength.txt
           3.  4.0K chrName.txt
           4.  4.0K chrStart.txt
           5.  3.0G Genome
           6.  4.0K genomeParameters.txt
           7.  1.2G Homo_sapiens.GRCh38.79.gtf
           8.  3.0G Homo_sapiens.GRCh38.dna.primary_assembly.fa
           9.  40K Log.out
          10.  23G SA
          11.  1.5G SAindex
  '''
  os.system(f'STAR --runThreadN {nthreads} --runMode genomeGenerate --genomeDir {idir} --genomeFastaFiles {idir}{fin}')
  print(f'STAR genome indexing completed successfully:\n\t{idir}{fin}')
