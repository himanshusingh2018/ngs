library(clusterProfiler)
library(enrichplot)
library(msigdbr)

dge <- read.csv('~/Downloads/DE_Significant_Results.csv')
dge <- dge[(dge$avg_log2FC>1 | dge$avg_log2FC < -1) & dge$p_val_adj < 0.05, ]

gene_list <- dge$avg_log2FC
names(gene_list) <- dge$X

senescence_genes <- msigdbr(species = "Homo sapiens",
                            category = "C2",
                            subcategory = "C2") %>%
  filter(gs_name == "FRIDMAN_SENESCENCE_UP") %>%
  dplyr::select(gene_symbol)

gse <- GSEA(geneList = gene_list,
            TERM2GENE = senescence_genes,
            valueCutoff = 0.05)

senescence_genes <- msigdbr(species = "Homo sapiens",collection="C2",subcollection = "CGP") %>%
  dplyr::filter(gs_name == c("FRIDMAN_SENESCENCE_UP",
                      "FRIDMAN_SENESCENCE_DN",
                      "TANG_SENESCENCE_TP53_TARGETS_UP",
                      "TANG_SENESCENCE_TP53_TARGETS_DN",
                      "KAMMINGA_SENESCENCE") %>%
  dplyr::select(gs_name, gene_symbol)

senescence_gene_sets <- msigdbr(species = "Homo sapiens", 
                                category = "C2", 
                                subcategory = "CGP") %>% 
  filter(gs_name %in% c("KAMMINGA_SENESCENCE", 
                        "TANG_SENESCENCE_TP53_TARGETS_UP", 
                        "TANG_SENESCENCE_TP53_TARGETS_DN", 
                        "FRIDMAN_SENESCENCE_DN", 
                        "FRIDMAN_SENESCENCE_UP")) %>% 
  dplyr::select(gs_name, gene_symbol)
  
  
gse <- GSEA(geneList = gene_list,
            TERM2GENE = senescence_genes,
            valueCutoff = 0.1)
