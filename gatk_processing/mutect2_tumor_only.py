import os
def mutect2_nofilter(idir, odir, genome_fa, sample, jvm=5, nthreads=4):
    '''
    GATK Mutect2 Variant Calling Tumor Only Mode [No Filter]
    https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

    RUN: 
        mutect2_nofilter(idir='/path/idir/', odir='/path/odir/', genome_fa = '/path/refgenome.fa', sample = 'samp_name', jvm=5, nthreads=4)
            idir     : '/path/idir/', 
            odir     : '/path/odir/', 
            genome_fa: '/path/refgenome.fa', 
            sample   : 'samp_name' 
            jvm      : java heapsize memory
            nthreads : number of threads
    '''
    os.makedirs(odir, exist_ok=True)#create odir if not exists
    print(f'gatk Mutect2 --java-options -Xmx{jvm}g -R {genome_fa} -I {idir}{sample}.mkdup.bam -O {odir}{sample}.vcf.gz --native-pair-hmm-threads {nthreads}')

def mutect2_filter(idir, odir, genome_fa, afonlygnomadvcf, hg38pon1000Gvcf, sample, jvm, nthreads):
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
        mutect2_filter(idir='/path/idir/', odir='/path/odir/', genome_fa = '/path/refgenome.fa', sample = 'samp_name', jvm=5, nthreads=4)
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
    os.system(f'gatk Mutect2 --java-options -Xmx{jvm}g -R {genome_fa} -I {idir}{sample}.mkdup.bam -O {odir}{sample}.vcf.gz --germline-resource {afonlygnomadvcf} --panel-of-normals {hg38pon1000Gvcf} --native-pair-hmm-threads {nthreads}')
    print('Successfully completed:\n\t{odir}{sample}.vcf.gz')

def filter_variant(idir, odir, genome_fa, sample):
    '''
    Filter Variants Using Mutect2

    filter_variant(idir='/path/idir/', odir='/path/odir/', genome_fa='/path/genome.fa', sample='sample_name')
        Input
            idir     : /path/inputdir/vcf/
            odir     : /path/outdir/vcf/
            genome_fa: /path/genome.fa
            sample   : sample.name
        Output
            /path/output/dir/var.filter.vcf.gz
    '''
    os.system('gatk FilterMutectCalls -V {idir}{sample}.vcf.gz -O {odir}{sample}.filter.vcf.gz -R {genome_fa}')
    print(f'Output: {odir}{sample}.filter.vcf.gz')
