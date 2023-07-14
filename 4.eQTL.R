###############
# Description #
###############
#* The purpose of this script is to examine regional effect on gene expression.
#* This script extracts eQTL profiles from the eQTL catalog. User can specify 
#* relevant  gene expression profile and from which tissue to examine.
#* 
#* Region of interest is dependent on locus while gene expression and tissue is
#* dependent on disease.
#* 
#* Genes to examine can be informed by using identified putative causal variants
#* identified in 3.statistical finemap script. It is also advised to use online 
#* databases such as vannoportal http://www.mulinlab.org/vportal/ to identify 
#* genes from putative causal variants
#* 

#* Steps
#* 1. Set up link with online eQTL catalog
#* 2. Define QTL region
#* 3. Extract eQTL for different tissue types
#* 4. Create Manhatten plots for each tissue
#* 5. Get top SNP

## Author: Dr. Owen Williams
##
## Date Created: 13-07-2023
##
## Email: owen.williams8@nhs.net

####################
# Install packages #
####################

require(tidyverse)
require(readr)
require(seqminer)
require(ggpubr)
require(GenomicRanges)
require(MungeSumstats)

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

########################
# 2. Define QTL region #
########################
#* Region of interest is on Chromsome 17 between BP 39500000, 40250000 however
#* to examine for cis regulatory effects we will expand region to 37000000 to 
#* 43000000.
#* 
#* All eQTL profiles use the GRCh38 whereas the original summary statistics file
#* uses GRch37. Therefore to be comparable to eQTL data will need to be 
#* converted from GRCh38 to GRCh37.

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

# Define region of interest in GRCh38 build
region = "17:32332223-32782223"

#identify tissues in different studies
GTEx_studies = dplyr::filter(tabix_paths, study_label == "GTEx")
Schmiedel_studies = dplyr::filter(tabix_paths, study_label == "Schmiedel_2018")
Alasoo_studies = dplyr::filter(tabix_paths, study_label == "Alasoo_2018")

##############################################
# 3. Extract eQTL for different tissue types #
##############################################

# 3.1 Create function -----------------------------------------------------
#* Create a function that can extracts eQTL data from eQTL catalog which uses 
#* GRCh38 build and converts to GRCh37.
#* 
#* Inputs required are tissue, Study, gene ensemble code, gene name and region of interest

extract_eQTL_tissue_GRCh38toGRCh37 = function(Tissue, STUDY, GENEofInterestEnsemble, geneName, GRCh38Region){
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
  # plot
  plot1 = ggplot(summary_stats1, aes(x = BP, y = -log10(P))) +
    labs(title = paste('tissue:',Tissue, "Gene:", geneName)) +
    geom_point() +
    theme_bw() +
    ylim(0,max(-log(summary_stats1$P,10)))
  list = list(plot1, summary_stats1)
  return(list)
}



#############################################
# 4. Create Manhatten plots for each tissue #
#############################################

# 3.2 Run function on specific tissues ---------------------------------------
# to be used in figure


# CD4+ --------------------------------------------------------------------

CD4_IKZF3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD4_T-cell_naive",
                                              STUDY = 'Schmiedel_2018',
                                              GENEofInterest = 'ENSG00000161405',
                                              geneName = 'IKZF3',
                                              GRCh38Region = "17:38843747-44922632")

CD4_GSDMB = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD4_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000073605',
                                               geneName = 'GSDMB',
                                               GRCh38Region = "17:38843747-44922632")

CD4_ORLD3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD4_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000172057',
                                               geneName = 'ORLD3',
                                               GRCh38Region = "17:38843747-44922632")

CD4_GSDMA = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD4_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000167914',
                                               geneName = 'GSDMA',
                                               GRCh38Region = "17:38843747-44922632")

CD4_PGAP3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD4_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000161395',
                                               geneName = 'PGAP3',
                                               GRCh38Region = "17:38843747-44922632")

CD4_MSL1 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD4_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000188895',
                                               geneName = 'MSL1',
                                               GRCh38Region = "17:38843747-44922632")

CD4_PNMT = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD4_T-cell_naive",
                                              STUDY = 'Schmiedel_2018',
                                              GENEofInterest = 'ENSG00000141744',
                                              geneName = 'PNMT',
                                              GRCh38Region = "17:38843747-44922632")


# CD8+ --------------------------------------------------------------------


CD8_IKZF3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD8_T-cell_naive",
                                              STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000161405',
                                              geneName = 'IKZF3',
                                               GRCh38Region = "17:38843747-44922632")

CD8_GSDMB = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD8_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000073605',
                                               geneName = 'GSDMB',
                                               GRCh38Region = "17:38843747-44922632")

CD8_ORLD3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD8_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000172057',
                                               geneName = 'ORLD3',
                                               GRCh38Region = "17:38843747-44922632")


CD8_GSDMA = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD8_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000167914',
                                               geneName = 'GSDMA',
                                               GRCh38Region = "17:38843747-44922632")

CD8_PGAP3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD8_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000161395',
                                               geneName = 'PGAP3',
                                               GRCh38Region = "17:38843747-44922632")

CD8_MSL1 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD8_T-cell_naive",
                                               STUDY = 'Schmiedel_2018',
                                               GENEofInterest = 'ENSG00000188895',
                                               geneName = 'MSL1',
                                               GRCh38Region = "17:38843747-44922632")

CD8_PNMT = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "CD8_T-cell_naive",
                                              STUDY = 'Schmiedel_2018',
                                              GENEofInterest = 'ENSG00000141744',
                                              geneName = 'PNMT',
                                              GRCh38Region = "17:38843747-44922632")

# Macrophage --------------------------------------------------------------


macrophage_IKZF3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "macrophage_naive",
                                               STUDY = 'Alasoo_2018',
                                               GENEofInterest = 'ENSG00000161405',
                                               geneName = 'IKZF3',
                                               GRCh38Region = "17:38843747-44922632")

macrophage_GSDMB = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "macrophage_naive",
                                               STUDY = 'Alasoo_2018',
                                               GENEofInterest = 'ENSG00000073605',
                                               geneName = 'GSDMB',
                                               GRCh38Region = "17:38843747-44922632")

macrophage_ORLD3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "macrophage_naive",
                                                      STUDY = 'Alasoo_2018',
                                                      GENEofInterest = 'ENSG00000172057',
                                                      geneName = 'ORLD3',
                                                      GRCh38Region = "17:38843747-44922632")

macrophage_GSDMA = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "macrophage_naive",
                                                      STUDY = 'Alasoo_2018',
                                                      GENEofInterest = 'ENSG00000167914',
                                                      geneName = 'GSDMA',
                                                      GRCh38Region = "17:38843747-44922632")

macrophage_PGAP3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "macrophage_naive",
                                                      STUDY = 'Alasoo_2018',
                                                      GENEofInterest = 'ENSG00000161395',
                                                      geneName = 'PGAP3',
                                                      GRCh38Region = "17:38843747-44922632")

macrophage_MSL1 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "macrophage_naive",
                                                      STUDY = 'Alasoo_2018',
                                                      GENEofInterest = 'ENSG00000188895',
                                                      geneName = 'MSL1',
                                                      GRCh38Region = "17:38843747-44922632")
macrophage_PNMT = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "macrophage_naive",
                                                     STUDY = 'Alasoo_2018',
                                                     GENEofInterest = 'ENSG00000141744',
                                                     geneName = 'PNMT',
                                                     GRCh38Region = "17:38843747-44922632")


# Lung --------------------------------------------------------------------



lung_IKZF3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "lung",
                                               STUDY = 'GTEx',
                                               GENEofInterest = 'ENSG00000161405',
                                               geneName = 'IKZF3',
                                               GRCh38Region = "17:38843747-44922632")

lung_GSDMB = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "lung",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000073605',
                                                geneName = 'GSDMB',
                                                GRCh38Region = "17:38843747-44922632")

lung_ORLD3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "lung",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000172057',
                                                geneName = 'ORLD3',
                                                GRCh38Region = "17:38843747-44922632")



lung_GSDMA = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "lung",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000167914',
                                                geneName = 'GSDMA',
                                                GRCh38Region = "17:38843747-44922632")

lung_PGAP3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "lung",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000161395',
                                                geneName = 'PGAP3',
                                                GRCh38Region = "17:38843747-44922632")

lung_MSL1 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "lung",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000188895',
                                                geneName = 'MSL1',
                                                GRCh38Region = "17:38843747-44922632")


lung_PNMT = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "lung",
                                               STUDY = 'GTEx',
                                               GENEofInterest = 'ENSG00000141744',
                                               geneName = 'PNMT',
                                               GRCh38Region = "17:38843747-44922632")

# esophagus mucosa --------------------------------------------------------



esophagus_mucosa_IKZF3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_mucosa",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000161405',
                                                geneName = 'IKZF3',
                                                GRCh38Region = "17:38843747-44922632")


esophagus_mucosa_GSDMB = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_mucosa",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000073605',
                                                geneName = 'GSDMB',
                                                GRCh38Region = "17:38843747-44922632")

esophagus_mucosa_ORLD3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_mucosa",
                                                            STUDY = 'GTEx',
                                                            GENEofInterest = 'ENSG00000172057',
                                                            geneName = 'ORLD3',
                                                            GRCh38Region = "17:38843747-44922632")



esophagus_mucosa_GSDMA = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_mucosa",
                                                            STUDY = 'GTEx',
                                                            GENEofInterest = 'ENSG00000167914',
                                                            geneName = 'GSDMA',
                                                            GRCh38Region = "17:38843747-44922632")


esophagus_mucosa_PGAP3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_mucosa",
                                                            STUDY = 'GTEx',
                                                            GENEofInterest = 'ENSG00000161395',
                                                            geneName = 'PGAP3',
                                                            GRCh38Region = "17:38843747-44922632")

esophagus_mucosa_MSL1 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_mucosa",
                                                            STUDY = 'GTEx',
                                                            GENEofInterest = 'ENSG00000188895',
                                                            geneName = 'MSL1',
                                                            GRCh38Region = "17:38843747-44922632")

esophagus_mucosa_PNMT = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_mucosa",
                                                           STUDY = 'GTEx',
                                                           GENEofInterest = 'ENSG00000141744',
                                                           geneName = 'PNMT',
                                                           GRCh38Region = "17:38843747-44922632")


# esophagus muscularis ----------------------------------------------------


esophagus_muscularis_IKZF3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_muscularis",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000161405',
                                                geneName = 'IKZF3',
                                                GRCh38Region = "17:38843747-44922632")


esophagus_muscularis_GSDMB = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_muscularis",
                                                STUDY = 'GTEx',
                                                GENEofInterest = 'ENSG00000073605',
                                                geneName = 'GSDMB',
                                                GRCh38Region = "17:38843747-44922632")

esophagus_muscularis_ORLD3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_muscularis",
                                                                STUDY = 'GTEx',
                                                                GENEofInterest = 'ENSG00000172057',
                                                                geneName = 'ORLD3',
                                                                GRCh38Region = "17:38843747-44922632")


esophagus_muscularis_GSDMA = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_muscularis",
                                                                STUDY = 'GTEx',
                                                                GENEofInterest = 'ENSG00000167914',
                                                                geneName = 'GSDMA',
                                                                GRCh38Region = "17:38843747-44922632")


esophagus_muscularis_PGAP3 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_muscularis",
                                                                STUDY = 'GTEx',
                                                                GENEofInterest = 'ENSG00000161395',
                                                                geneName = 'PGAP3',
                                                                GRCh38Region = "17:38843747-44922632")

esophagus_muscularis_MSL1 = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_muscularis",
                                                                STUDY = 'GTEx',
                                                                GENEofInterest = 'ENSG00000188895',
                                                                geneName = 'MSL1',
                                                                GRCh38Region = "17:38843747-44922632")

esophagus_muscularis_PNMT = extract_eQTL_tissue_GRCh38toGRCh37(Tissue =  "esophagus_muscularis",
                                                               STUDY = 'GTEx',
                                                               GENEofInterest = 'ENSG00000141744',
                                                               geneName = 'PNMT',
                                                               GRCh38Region = "17:38843747-44922632")



CD4_IKZF3[[1]]
CD4_GSDMB[[1]]
CD4_ORLD3[[1]]
CD4_GSDMA[[1]]
CD4_PGAP3[[1]]
CD4_MSL1[[1]]
CD4_PNMT[[1]]

CD8_IKZF3[[1]]
CD8_GSDMB[[1]]
CD8_ORLD3[[1]]
CD8_GSDMA[[1]]
CD8_PGAP3[[1]]
CD8_MSL1[[1]]
CD8_PNMT[[1]]

macrophage_IKZF3[[1]]
macrophage_GSDMB[[1]]
macrophage_ORLD3[[1]]
macrophage_GSDMA[[1]]
macrophage_PGAP3[[1]]
macrophage_MSL1[[1]]
macrophage_PNMT[[1]]

lung_IKZF3[[1]]
lung_GSDMB[[1]]
lung_ORLD3[[1]]
lung_GSDMA[[1]]
lung_PGAP3[[1]]
lung_MSL1[[1]]
lung_PNMT[[1]]

esophagus_mucosa_IKZF3[[1]]
esophagus_mucosa_GSDMB[[1]]
esophagus_mucosa_ORLD3[[1]]
esophagus_mucosa_GSDMA[[1]]
esophagus_mucosa_PGAP3[[1]]
esophagus_mucosa_MSL1[[1]]
esophagus_mucosa_PNMT[[1]]

esophagus_muscularis_IKZF3[[1]]
esophagus_muscularis_GSDMA[[1]]
esophagus_muscularis_GSDMB[[1]]
esophagus_muscularis_ORLD3[[1]]
esophagus_muscularis_PGAP3[[1]]
esophagus_muscularis_MSL1[[1]]
esophagus_muscularis_PNMT[[1]]


# Create plot -------------------------------------------------------------

ggarrange(CD4_IKZF3[[1]], CD4_GSDMB[[1]], CD4_ORLD3[[1]], CD4_GSDMA[[1]], CD4_PGAP3[[1]], CD4_MSL1[[1]], NULL,
          CD8_IKZF3[[1]],CD8_GSDMB[[1]],CD8_ORLD3[[1]],CD8_GSDMA[[1]],CD8_PGAP3[[1]],CD8_MSL1[[1]],CD8_PNMT[[1]],
          NULL,macrophage_GSDMB[[1]],macrophage_ORLD3[[1]],macrophage_GSDMA[[1]],macrophage_PGAP3[[1]],macrophage_MSL1[[1]],NULL,
          ncol = 7, nrow = 3)


ggarrange(lung_IKZF3[[1]],lung_GSDMB[[1]],lung_ORLD3[[1]],lung_GSDMA[[1]],lung_PGAP3[[1]],lung_MSL1[[1]],lung_PNMT[[1]],
          esophagus_mucosa_IKZF3[[1]],esophagus_mucosa_GSDMB[[1]],esophagus_mucosa_ORLD3[[1]],esophagus_mucosa_GSDMA[[1]],esophagus_mucosa_PGAP3[[1]],esophagus_mucosa_MSL1[[1]],esophagus_mucosa_PNMT[[1]],
          esophagus_muscularis_IKZF3[[1]],esophagus_muscularis_GSDMA[[1]],esophagus_muscularis_GSDMB[[1]],esophagus_muscularis_ORLD3[[1]],esophagus_muscularis_PGAP3[[1]],esophagus_muscularis_MSL1[[1]],esophagus_muscularis_PNMT[[1]],
          ncol = 7, nrow = 3)

#######################
# 5. Get leading SNPs #
#######################

#* Get LeadSNPs for each tissue. In this example I have only done this for lung,
#* however this should be repeated for all tissues
#* 

LungleadSNPs = lung_IKZF3[[2]] %>%
  dplyr::filter(P == min(P)) %>%
  mutate(Tissue = 'lung', geneName = 'IKZF3')

LungleadSNPs = rbind(LungleadSNPs,
                 lung_GSDMB[[2]] %>%
                   dplyr::filter(P == min(P)) %>%
                   mutate(Tissue = 'lung', geneName = 'GSDMB'),
                 lung_ORLD3[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue ='lung', geneName = 'ORLD3'),
                 lung_GSDMA[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue ='lung', geneName = 'GSDMA'),
                 lung_PGAP3[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue ='lung', geneName = 'PGAP3'),
                 lung_MSL1[[2]] %>%
                   dplyr::filter(P == min(P)) %>%
                   mutate(Tissue ='lung', geneName = 'MSL1'),
                 lung_PNMT[[2]] %>%
                   dplyr::filter(P == min(P))%>%
                   mutate(Tissue ='lung', geneName = 'PNT'))


# remove any lead SNPs where below the -log10(p) = 5 threshold

LungleadSNPs = LungleadSNPs %>%
  filter(-log10(P) > 5)

# Export table

write_csv(LungleadSNPs, 'Results/eQTLLeadSNPs.csv')
