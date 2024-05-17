def read_vcf_as_dataframe(fvcf):
    '''
    Read vcf
        Comment all rows starting as '#'
        Column Names: starting with '#CHROM'
    '''
    col = next(l for l in open(fvcf) if l.startswith('#CHROM')).strip().split('\t')
    vcf = pd.read_csv(fvcf, sep="\t",header=None, comment='#', names=col)
    return(vcf)

read_vcf_as_dataframe('../wes_analysis/PID_14/sv/delly/pid14_3.filter.vcf')
