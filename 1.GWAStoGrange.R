###############
# Description #
###############
#* Standardise imported GWAS file which is compatible with GenomicRanges

## Author: Dr. Owen Williams
##
## Date Created: 12-07-2023
##
## Email: owen.williams8@nhs.net


####################
# Install packages #
####################
require(tidyverse)
require(MungeSumstats)
require(GenomicRanges)
require(sqldf)

#################
# Set directory #
#################

here::here()

##########################################################
#  Download GWAS Data sets and convert to Genomic Ranges #
##########################################################
# Summary Statistics previously downloaded 

# 1. Read in file with SQL query
#* add sep = /t to covert to tsv

GWAS = sqldf::read.csv.sql(file = 'Raw data/29273806-GCST006862-EFO_0000270.h.tsv',
                           sep = '\t')


colnames(GWAS)

# wrangle data
GWAS = GWAS %>%
  dplyr::select(hm_chrom, hm_pos, hm_rsid, hm_beta, p_value, standard_error)

# Change col names
colnames(GWAS) = c('CHR','BP', 'SNP', 'BETA', 'P', 'SE')

# Remove any SNPs with unspecified Chromosome and arrange by ascending order
GWAS = GWAS  %>%
  filter(CHR != 0) %>%
  arrange(CHR, BP) %>%
  filter(duplicated(SNP) == FALSE)

# modify Chromosome column to add 'chr'
GWAS$CHR = paste0('chr', GWAS$CHR)

# Create GenomicRanges file
GWAS_genomeR = makeGRangesFromDataFrame(GWAS,
                                           seqnames.field=c("CHR"),
                                           start.field="BP",
                                           end.field=c("BP"))

# Add p-values, SE and BETA
GWAS_genomeR$P = GWAS$P
GWAS_genomeR$SE = GWAS$SE
GWAS_genomeR$BETA = GWAS$BETA
GWAS_genomeR$SNP = GWAS$SNP
# rsIDs
names(GWAS_genomeR) = GWAS$SNP

# view genomic range file
GWAS_genomeR

# Remove dataframe to save memory
rm(GWAS)

# write Grange to file
write.table(as.data.frame(GWAS_genomeR),
            file="Grange/Granges_GWASGRCh37.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")
