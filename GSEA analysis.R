###############
# Description #
###############
#* The purpose of this script is to perform a gene set enrichment analysis 
#* (GSEA) on co-expressed genes. This includes:
#* 
#* 1. GO enrichment GSEA
#* 2. KEGG enrichment GSEA
#* 
#*  Outputs of this GSEA may provide further insights into druggable targets.
#* 
#* In this example we are using co-expressed genes to previously identified
#* putative causal gene to asthma: ORMDL3. The list of co-expressed genes was
#* downloaded from the https://coxpresdb.jp/gene_coexpression/ website.
#* 

## Author: Dr. Owen Williams
##
## Date Created: 14-07-2023
##
## Email: owen.williams8@nhs.net

####################
# Install packages #
####################

# Load the biomaRt library
require(clusterProfiler)
require(org.Hs.eg.db)
require(tidyverse)

#########################
# Conduct GSEA analysis #
#########################

# upload co-expressed gene list
coexpression_ORMDL3 = read_csv(file = 'Raw data/coexpression.csv')

# Perform GO enrichment analsysis
GO_GSEA_result <- enrichGO(
  gene         = coexpression_ORMDL3$Gene,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",  # biological process
  pAdjustMethod = "fdr",  # false detection rate
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

# See results
as.data.frame(GO_GSEA_result)

# Perform KEGG enrichment analysis
kegg_GSEA_result <- enrichKEGG(
  gene         = coexpression_ORMDL3$`Entrez Gene ID`,
  organism     = 'hsa', 
  keyType      = 'kegg',
  pvalueCutoff = 0.05,  
  pAdjustMethod = "fdr")

# See results
as.data.frame(kegg_GSEA_result)

