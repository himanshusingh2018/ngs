import os, pandas as pd, numpy as np,sys

def treeomics_input_by_impact_gene_panel():
    '''   
    COMMAND: python treeomics_input_by_impact_gene_panel.py idir odir
        idir : /path/input/dir: takes input files: vcf/coverage.txt; vcf/mut_read_table.txt
        odir : /path/output/dir: output will be generated: vcf/coverage.impact.gene.panel.txt;vcf/mut_read_table.impact.gene.panel.txt 
    e.g.
    python treeomics_input.filter.impact.panel.py PID_113 PID_113 impact_gene_panel/IMPACT341-Gene-List_20140101.xlsx IMPACT341
    '''
    idir = sys.argv[1]; odir = sys.argv[2]; gene_panel_excel_file = sys.argv[3]; sheet = sys.argv[4]

    gp = pd.read_excel(gene_panel_excel_file, sheet_name=sheet)#read gene panel
    cov = pd.read_csv(f'{idir}/vcf/coverage.txt',sep="\t",header=0)#read coverage.txt
    mut = pd.read_csv(f'{idir}/vcf/mut_read_table.txt',sep="\t",header=0)#read mut_read_table.txt
    #writing files
    cov[cov.Gene.isin(gp.Gene_Symbol)].to_csv(f'{odir}/vcf/coverage.{sheet}.txt', sep="\t", index=False)
    mut[mut.Gene.isin(gp.Gene_Symbol)].to_csv(f'{odir}/vcf/mut_read_table.{sheet}.txt', sep="\t", index=False)
    
    print(f"Successfully generated:\n\t1: {odir}/vcf/coverage.{sheet}.txt\n\t2: {odir}/coverage.{sheet}.txt")

treeomics_input_by_impact_gene_panel()
#python treeomics_input.filter.impact.panel.py PID_113 PID_113 impact_gene_panel/IMPACT341-Gene-List_20140101.xlsx IMPACT341