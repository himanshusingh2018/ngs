import os
def picard_bam_clean(idir, odir, sample):
    '''
    PICARD CLEANING SORT BAM FILE
    
    Reason: https://davetang.org/wiki/tiki-index.php?page=SAM#Clipped_alignment
    
            One query sequence may be aligned to multiple places on the reference genome, either with or without overlaps.
            Hard clipping: Avoid presenting the full query sequence multiple times for non-overlapping hits
            Soft clipping: these sequence are not part of the alignment. 
            Therefore, soft clipping shall be removed. 
    RUN:
        picard_bam_clean(idir = '/path/to/bam/file/', odir='/path/to/output/dir/', sample=sample_name)

    OUTPUT: /path/to/output/dir/sample.clean.bam

    Picard CleanSam Arguments:
        cmd: picard CleanSam --INPUT /bam/file/path/sample.sort.bam --OUTPUT /output/path/sample.clean.bam

                CleanSam: Cleans the BAM, soft-clipping beyond-end-of-reference alignments and setting MAPQ (mapping quality) to 0 for unmapped reads 
                            soft-clipping: bases in 5' and 3' of the read are NOT part of the alignment
                            hard-clipping: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been removed from the read sequence in the BAM file.
                INPUT   : Input, sort bam file
                OUTPUT  : Output, clean bam file

    '''
    os.makedirs(odir, exist_ok=True)#creat odir if not exists
    os.system(f"picard CleanSam --INPUT {idir}{sample}.sort.bam --OUTPUT {odir}{sample}.clean.bam")
    print(f"{odir}{sample}.clean.bam is generated successfully...")

def picard_markduplicate(idir, odir, sample):
    '''
    Mark duplicate reads.
    Setting MAPQ to 0 for unmapped reads

    Reason: duplicated mark reads are ignored by the post-processing. The duplicated reads should not removed, as they may impact the coverage.
    
    RUN: 
        picard_markduplicate(idir='/path/to/bam/', odir='/path/to/bam/', sample='sample_name.clean.bam')
    
    OUTPUT: Two output files:
                1. {/path/to/output/}{sample}.sort.mkdup.bam
                2. {/path/to/output/}{sample}.sort.dup.metrics.txt

    Picard Markduplicate Arguments:
        --INPUT                 : input SAM, BAM or CRAM files to analyze. Must be coordinate sorted.
        --OUTPUT                : output file to write marked records
        --METRICS_FILE          : File to write duplication metrics to  Required.  
        --VALIDATION_STRINGENCY : {STRICT, LENIENT, SILENT}. LENIENT Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT.
        --CREATE_INDEX          : {true|false}, Whether to create an index when writing VCF or coordinate sorted BAM output.  Default value: false. Possible values: {true, false} 
    '''
    os.makedirs(odir, exist_ok=True)#creat odir if not exists
    os.system(f"picard MarkDuplicates --INPUT {idir}{sample}.clean.bam --OUTPUT {odir}{sample}.sort.mkdup.bam --METRICS_FILE {odir}{sample}.sort.dup.metrics.txt --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true")
    print(f"Successfully generated files: \n\t1: {odir}{sample}.sort.mkdup.bam \n\t2: {odir}{sample}.sort.dup.metrics.txt")
