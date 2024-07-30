import pandas as pd, gzip, sys

def get_col_names(file_path):
    """
    Extract column names from a VCF file.
    
    Args:
    - file_path (str): Path to the VCF file.

    Returns:
    - list: List of column names.
    """
    with (gzip.open(file_path, 'rt', encoding='ISO-8859-1') if file_path.endswith('.gz') else open(file_path, 'r', encoding='ISO-8859-1')) as f:
        # Find and return the column names from the header line starting with '#CHROM'
        return next(l.strip().split('\t') for l in f if l.startswith('#CHROM'))

def process_funcotator_vcf(file_path):
    """
    Process a Funcotator VCF file.
    
    Args:
    - file_path (str): Path to the VCF file.

    Returns:
    - pd.DataFrame: Processed DataFrame with specific columns and filters applied.
    """
    col = get_col_names(file_path)
    
    # Read VCF file
    vcf = pd.read_csv(
        gzip.open(file_path, 'rt', encoding='ISO-8859-1') if file_path.endswith('.gz') else open(file_path, 'r', encoding='ISO-8859-1'),
        sep="\t", header=None, comment='#', names=col
    )

    # Filter rows where FILTER column is 'PASS'
    vcf = vcf[vcf.FILTER == 'PASS']

    # Process INFO column and expand 'New' into separate columns
    vcf['New'] = vcf['INFO'].apply(lambda x: '|'.join(x.split(';')[4].split('|')[:10]) if len(x.split(';')) > 2 else None)
    vcf['New'] = vcf['New'].apply(lambda x: x.replace('FUNCOTATION=[', '') if pd.notna(x) else x)
    vcf[['Symbol', 'HgVer', 'Chr', 'Start', 'End', 'Feature', 'Empty', 'Mut', 'Ref', 'Alt']] = vcf['New'].str.split('|', expand=True)

    # Extract relevant columns and split additional information
    req_col = ['Symbol', 'Chr', 'Start', 'End', 'Feature', 'Mut', 'REF', 'ALT'] + col[9:]
    for x in req_col[8:]:
        vcf[x] = vcf[x].apply(lambda cell: cell.split(':')[0] if pd.notna(cell) else None)

    # Filter out specific features
    exclude_feature = ['INTRON', 'IGR', 'FIVE_PRIME_UTR', 'THREE_PRIME_UTR', 'LINCRNA']
    vcf_filtered = vcf[~vcf['Feature'].isin(exclude_feature)][req_col]

    return vcf_filtered

# Example usage
#fvcf = '/Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/PID_1/vcf/pid1.filter.funcotator.vcf'
fvcf = sys.argv[1]
filtered_vcf = process_funcotator_vcf(fvcf)
print(filtered_vcf)

#run the script
#python oncoprint_data_funcotator.py <funcotator_annot.vcf>
#python oncoprint_data_funcotator.py /Volumes/lilac_data_ziv/transciptome/paired_pnet/wes_analysis/PID_1/vcf/pid1.filter.funcotator.vcf