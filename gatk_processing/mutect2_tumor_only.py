import os

def tumoronly_nofilter(idir, odir, genome_fa, sample, jvm=5, nthreads=4):
    '''
    GATK Mutect2 Variant Calling Tumor Only Mode [No Filter]
    https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

    RUN: 
        tumoronly_nofilter(idir='/path/idir/', odir='/path/odir/', genome_fa = '/path/refgenome.fa', sample = 'samp_name', jvm=5, nthreads=4)
            idir     : '/path/idir/', 
            odir     : '/path/odir/', 
            genome_fa: '/path/refgenome.fa', 
            sample   : 'samp_name' 
            jvm      : java heapsize memory
            nthreads : number of threads
    '''
    os.makedirs(odir, exist_ok=True)#create odir if not exists
    print(f'gatk Mutect2 --java-options -Xmx{jvm}g -R {genome_fa} -I {idir}{sample}.bam -O {odir}{sample}.vcf.gz --native-pair-hmm-threads {nthreads}')

def tumoronly_filter(idir, odir, genome_fa, afonlygnomadvcf, hg38pon100gvcf, sample, jvm=5, nthreads=4):
    '''
    GATK Mutect2 Variant Calling Tumor Only Mode [FILTER]
    https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
    
        Filter by: (https://github.com/gwcbi/genomics_prostate_cancer/blob/master/refs/protocol.sh)
            1. germline-resource:
                Download gnomad
                    wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
                    wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
 
            2. panel of normals:
                wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
                wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi

    RUN: 
        tumoronly_nofilter(idir='/path/idir/', odir='/path/odir/', genome_fa = '/path/refgenome.fa', sample = 'samp_name', jvm=5, nthreads=4)
                idir            : '/path/idir/', 
                odir            : '/path/odir/', 
                genome_fa       : '/path/refgenome.fa', 
                sample          : 'samp_name' 
                jvm             : java heapsize memory
                nthreads        : number of threads
                afonlygnomadvcf : af-only-gnomadvcf.vcf.gz  
                hg38pon1000gvcf  : 1000g_pon.hg38.vcf.gz

    '''
    os.makedirs(odir, exist_ok=True)#create odir if not exists
    print(f'gatk Mutect2 --java-options -Xmx{jvm}g -R {genome_fa} -I {idir}{sample}.bam -O {odir}{sample}.vcf.gz --germline-resource {afonlygnomadvcf} --panel-of-normals {hg38pon1000gvcf} -O {odir}{sample}.vcf.gz --native-pair-hmm-threads {nthreads}')


