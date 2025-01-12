import pandas as pd
from pyensembl import EnsemblRelease 
'''
python version: 3.9
numpy: 1.29
'''
# Load the genome annotation for human (GRCh38) 
data = EnsemblRelease(104, species='homo_sapiens') 
# Download and index the data (if not already done) 
data.download() 
data.index() 
# Function to get gene name and exon ID by coordinate 
def get_gene_and_exon_by_coordinate(chrom, pos):
    # Get genes at the locus 
    genes = data.genes_at_locus(contig=chrom, position=pos) 
    gene_names = [gene.name for gene in genes] 
    # Get exons at the locus 
    exons = data.exons_at_locus(contig=chrom, position=pos) 
    exon_ids = [exon.id for exon in exons] 
    return gene_names, exon_ids

t1 = pd.read_csv('~/Desktop/test1.txt',sep="\t",header=0)

for idx, row in t1.iterrows():
    chr = str(row.Chromosome)
    pos = row.Position
    res = get_gene_and_exon_by_coordinate(chr,pos)
    print(chr,pos,res)
