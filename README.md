# ngs
NGS Data Analysis

# used python libraries
    import os, sys, glob
    import pandas as pd
    import functools
    from functools import reduce


# project name
    cnv_analysis
        'Need to add modules'

    gatk_processing
        vcfannotation: VCF Annotation using Funcotator
        mutect2_tumor_only
            tumoronly_nofilter : Mutect2 Somatic Variant Calling Tumor Only Mode [No Filter]
            tumoronly_filter   : Mutect2 Somatic Variant Calling Tumor Only Mode [Filter: af-only-gnomad.vcf.gz and 1000g_pon.hg38.vcf.gz]
            filter variant     : Mutect2 Filter Variants [PASS]

    ngs_analysis
        genome_index
            bowtie2_gi: Bowtie2 Genome Indexing
            bwa_gi    : BWA Genome Indexing
            star_gi   : STAR Genome Indexing 

        read_align
            bowtie2_fastq_aln: bowtie2 FastQ Alignment (Not for variant calling)
            bwa_fastq_aln    : BWA MEM FastQ Alignment (preferred for Variant)

        star_rnaseq_fastq_2pass_aln: STAR 2 Pass FastQ Alignment    


        htseq_count
            htseq_exon_count     : HTSEQ Exon Count
            merge_all_htseqcount : Merege All htseqcount in One

        picard_process
            picard_bam_clean     : Picard Cleaning Sorted BAM file for soft clipping
            picard_markduplicate : Mark duplicates in bam, Unmapped reads MAPQ set 0

    treeomics_analysis
        treeomics_input
        
    misc
        extract_data_between_two_pattern
        gtf_process
        sample_vs_gene_mut_table