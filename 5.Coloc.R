###############
# Description #
###############
#* 
#* This script is used as a quality control check to determine whether the eQTL
#* and GWAS variant share same causal variant.
#* 
#* This will perform co-localisation analysis using lung eQTL profiles with the
#* genes where a significant eQTL which was performed in the eQTL.R script with
#* results saved eQTLLeadSNPs.csv
#* 
#* This includes eQTL data examines gene expresssion in:
#*  1. GSDMB (ENSG00000073605)
#*  2. ORLD3 (ENSG00000172057)
#*  3. GSDMA (ENSG00000167914)
#*  4. PGAP3 (ENSG00000161395)
#*  5. MSL1 (ENSG00000188895)
#* 
#* Steps
#* 1. Set up link with online eQTL catalog
#* 2. get eQTL data
#* 3. Convert to GRCH37 build
#* 4. Get GWAS study 
#* 5. Get LD matrix
#* 6. Compile datasets for coloc 
#* 7. Coloc for Single Causal variant
#* 8. Repeat for other GWAS studies

#* WARNING: These functions depend on LDlinkR API and can fail if server is
#* overloaded. If this happens repeat function until it works
#* 
#* WARNING: As this is working with large data some of these functions may take 
#* a long time populate.

## Author: Dr. Owen Williams
##
## Date Created: 13-07-2023
##
## Email: owen.williams8@nhs.net


####################
# Install packages #
####################

require(tidyverse)
library(readr)
library(coloc)
library(seqminer)
require(ggpubr)
require(LDlinkR)
require(MungeSumstats)
library(biomaRt)

###########################################
# 1. Set up link with online eQTL catalog #
###########################################

# Set link to online eQTL catalog
tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Get list of eQTL data from GTex study 
imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Create function for extracting information from study using Tabix
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  colnames(fetch_table) = column_names
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::select(-rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
}

###################
# 3. Extract eQTL #
###################


# define region of interest in GRCh37 build
sumstats_dt = data.frame(SNP = 1:2,
                         CHR = 17,
                         BP = c(37000000, 43000000))

# convert to GRCh38 build
sumstats = liftover(sumstats_dt=sumstats_dt, 
                    ref_genome = "hg19",
                    convert_ref_genome="hg38")

# Get new range in GRCh38 build
min(sumstats$BP)
max(sumstats$BP)


# Create function
eQTL_GRCh38toGRCh37_forColoc = function(Tissue, STUDY, GENEofInterestEnsemble, geneName, GRCh38Region){
  studies_df = dplyr::filter(tabix_paths, study_label == STUDY, sample_group == Tissue)
  #Extract column names from first file
  column_names = colnames(readr::read_tsv(studies_df$ftp_path, n_max = 1))
  #Import summary statistics - selected gene = HLA-DRB1
  summary_stats1 = import_eQTLCatalogue(studies_df$ftp_path, GRCh38Region, selected_gene_id = GENEofInterestEnsemble, column_names) 
  
  # Convert column names to standardised format
  summary_stats1 = summary_stats1 %>%
    dplyr::rename(CHR = chromosome) %>%
    dplyr::rename(SNP = variant) %>%
    dplyr::rename(BP = position) %>%
    dplyr::rename(P = pvalue)
  
  # Convert to GRCh38 to GRCh37 using mungestats lift over
  summary_stats1 = MungeSumstats::liftover(sumstats_dt = summary_stats1,
                                           ref_genome = "hg38",
                                           convert_ref_genome = "hg19")
  return(summary_stats1)
}

#################
# Get eQTL data #
#################

lung_GSDMB = eQTL_GRCh38toGRCh37_forColoc(Tissue =  "lung",
                                           STUDY = 'GTEx',
                                           GENEofInterest = 'ENSG00000073605',
                                           geneName = 'GSDMB',
                                           GRCh38Region = "17:38843747-44922632")
##################
# Get GWAS study #
##################

GWASlocus = read_tsv('Grange/Granges_GWASGRCh37.tsv') %>%
  filter(seqnames == 'chr17',
         between(start, 37000000, 43000000)) %>%
  select(seqnames, start, P, SE, BETA, SNP) %>%
  rename(CHR = seqnames,
         BP = start)

###################################################
# Filter both studies so they are below 1000 SNPs #
###################################################
#* Make sure both studies are below the 1000 SNP threshold so that SNP matrix 
#* can be obtained from LDLinkR website

GWASlocus_filter = GWASlocus %>%
  filter(-log10(P) > 1)

plot(GWASlocus_filter$BP,-log10(GWASlocus_filter$P))


lung_GSDMB_filter = lung_GSDMB %>%
  filter(-log10(P) > 3) %>%
  distinct(BP, .keep_all = TRUE)

plot(lung_GSDMB_filter$BP,-log10(lung_GSDMB_filter$P))


#################
# Get LD matrix #
#################

#LDlink token = b68c33126cbf
token = 'b68c33126cbf'

# Get LD matrix for GWAS data --------------------------------------------------

LD_GWAS = LDlinkR::LDmatrix(snps = GWASlocus_filter$SNP, pop = "EUR", r2d = "r2",
                       token = token, file = TRUE)



# Get LD Matrix for eQTL data ---------------------------------------------
#* Because eQTL does not come with rsid information which is required for LD link
#* we will get this from bioMart
#* 
#* 

# get genomic co-ordinates of eQTL variants
df = data.frame(CHR = as.character(lung_GSDMB_filter$CHR),
                chr_start = as.character(lung_GSDMB_filter$BP),
                chr_end = as.character(lung_GSDMB_filter$BP))

position <- apply(df, 1, paste, collapse = ":")
position

# Define the SNP mart
snp_mart <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Get RSID information using biomart
snp_data = getBM(attributes = c('refsnp_id', 'chrom_start'), 
                 filters = 'chromosomal_region', 
                 values = position, 
                 mart = snp_mart)

# tidy data
snp_data = snp_data %>%
  dplyr::rename(BP = chrom_start)

# left join RSID information to eQTL data
lung_GSDMB_filter_RSID = dplyr::left_join(lung_GSDMB_filter, snp_data, by = 'BP')

# remove any non RSids and duplicates
lung_GSDMB_filter_RSID = lung_GSDMB_filter_RSID %>%
  filter(!str_starts(refsnp_id, "rs")) %>%
  distinct(BP, .keep_all = TRUE)

# Get LD
LD_GSDMB = LDlinkR::LDmatrix(snps = lung_GSDMB_filter_RSID$refsnp_id, pop = "EUR", r2d = "r2",
                            token = token, file = TRUE)


# tidy GWAS LD matrix -----------------------------------------------------

# filter SNPS in GWAS which were not found in LD matrix
GWASlocus_filter = GWASlocus_filter %>%
  filter(SNP %in% LD_GWAS$RS_number)

# remove SNP names as column
LD_GWAS.x = LD_GWAS[-1]
# rownames = colnmaes
rownames(LD_GWAS.x) = GWASlocus_filter$BP
colnames(LD_GWAS.x) = rownames(LD_GWAS.x)
# Convert to matrix and remove NAs
LD_GWAS.x = data.matrix(LD_GWAS.x %>%
                     mutate(across(everything(), ~ replace_na(.x, 0))))

# tidy eQTL LD matrix ----------------------------------------------------------

# filter SNPS in eQTL which were not found in LD matrix
lung_GSDMB_filter_RSID = lung_GSDMB_filter_RSID %>%
  filter(refsnp_id %in% LD_GSDMB$RS_number)

# repeat LD extraction with new SNP list
LD_GSDMB = LDlinkR::LDmatrix(snps = lung_GSDMB_filter_RSID$refsnp_id, pop = "EUR", r2d = "r2",
                             token = token, file = TRUE)

# filter SNPS in LD which were not found in eQTL matrix
LD_GSDMB = LD_GSDMB %>%
  filter(RS_number %in% lung_GSDMB_filter_RSID$refsnp_id)

# remove SNP names as column
LD_GSDMB.x = LD_GSDMB[-1]
# rownames = colnmaes
rownames(LD_GSDMB.x) = lung_GSDMB_filter_RSID$BP
colnames(LD_GSDMB.x) = rownames(LD_GSDMB.x)
# Convert to matrix and remove NAs
LD_GSDMB.x[is.na(LD_GSDMB.x)] = 0

LD_GSDMB.x = data.matrix(LD_GSDMB.x)

######################################
# Coloc for Multiple causal variants #
######################################
#* Perform colocalisation analysis using the SuSiE regression framework that allows
#* multiple causal variants 

# Run Susie on both datasets
zscoreGWAS = GWASlocus_filter$BETA/GWASlocus_filter$SE
S1 = susieR::susie_rss(zscoreGWAS, data.matrix(LD_GWAS.x), n =1000)
summary(S1)


zscoreeQTL = lung_GSDMB_filter_RSID$beta/lung_GSDMB_filter_RSID$se
S2 = susieR::susie_rss(zscoreeQTL, data.matrix(LD_GSDMB.x), n =1000)
summary(S2)


# Perform colocalisation
x =
  
  if(requireNamespace("susieR",quietly=TRUE)) {
    susie.res=coloc.susie(S1,S2)
    print(susie.res$summary)
  }

write_csv(x, '/Users/owen/Desktop/ColocVariants.csv')

################################################################################
################## Repeat for different gene profiles ##########################
################################################################################