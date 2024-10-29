from pyensembl import EnsemblRelease
data = EnsemblRelease(75)# Choose the Ensembl release (e.g., 75 for GRCh37)

def find_gene(chromosome, position):#annotate gene name by coordinate
    try:
        genes = data.genes_at_locus(contig=str(chromosome), position=int(position))
        if genes:
            return genes[0].gene_name
        else:
            return "No gene found"
    except Exception as e:
        return f"Error: {str(e)}"

# Assuming df is your DataFrame with 'chromosome' and 'position' columns
voi_df['GENE'] = voi_df.apply(lambda row: find_gene(row.CHR[3:], row.POS), axis=1)

