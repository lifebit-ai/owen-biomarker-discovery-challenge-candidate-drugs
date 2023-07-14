###############
# Description #
###############
#* Fine-mapping causal variants using SuSiE finemapping package, from summary
#* statistics data, including importing linkage reference panel from LDlink
#* 
#* Steps
#* 1. Setup LD link to extract LD matrix
#* 2. Import and filter GWAS files
#* 3. Create Function to automate LD download and SuSiE finemapping
#* 4. Fine map GWAS files
#* 5. Create Function to plot LD heatmap
#* 6. Create LD map
#* 7. Summarise data

#* 6. Produce LDheat maps
#* 7. Concatenate plots

#* WARNING: These functions depend on LDlinkR API and can fail if server is
#* overloaded. If this happens repeat function until it works
#* 
#* WARNING: As this is working with large data some of these functions may take 
#* a long time populate.

## Author: Dr. Owen Williams
##
## Date Created: 12-07-2023
##
## Email: owen.williams8@nhs.net

####################
# Install packages #
####################
require(data.table)
require(tidyverse)
require(susieR)
require(LDlinkR)
require(gaston)

####################
# 1. Set up LDlink #
####################

#LDlink token = b68c33126cbf
token = 'b68c33126cbf'


##################################
# 1.Read in Gr_ranges if created #
##################################

# Read in file  -----------------------------------------------

Locus = read_tsv('Grange/Locus.tsv') %>%
  dplyr::select(seqnames, start, SNP, BETA, SE, pval)


# Standarise column names
colnames(Locus) = c('CHR', 'BP', 'SNP', 'BETA', 'SE', 'P')

# remove any SNPs without rs number
LocusSNP_list = Locus %>%
  dplyr::filter(str_detect(Locus$SNP, 'rs'))

# Plot region containing SNPs to be fine-mapped
plot(-log10(LocusSNP_list$P)~LocusSNP_list$BP)

############################################
# 3. Develop function to automate Finemapping #
############################################
#* Create a function which automates SuSiE fine mapping 
#* 
#* Using SNPs of interest this function will download suitable reference LD panel
#* from LDLink using specified population, calculate Z-score for each variant and 
#* perform SuSiE fine mapping to identify credible sets with 95% confidence of 
#* non zero effect.
#* 
#* Inputs:
#* 1) Loci  
#* 2) LD population used in summary statistic
#*        

#* Outputs include:
#* i) SuSiE_rss fit,   
#* ii) highlighted SNPs plot,
#* iii) table of credible sets
#* iv) info of credible sets

get_cc = function(Locus, LDLinkPopulation){
  #* Import LD matrix using LD link
  LD = LDlinkR::LDmatrix(snps = Locus$SNP, pop = LDLinkPopulation, r2d = "r2",
                         token = token, file = TRUE)
  
  # filter SNP list from inputted GWAS file and remove any missing variants 
  SNP_list = Locus %>%
    filter(SNP %in% as.vector(LD[1]$RS_number))
  
  # Convert LD into a matrix and convert NAs to 0 to make SuSiE compatabilble
  LDmatrix = data.matrix(LD[-1] %>%
                           mutate(across(everything(), ~ replace_na(.x, 0))))
  
  # specify how many samples included
  n=nrow(LD)
  
  # calculate z-score
  zscore = SNP_list$BETA/SNP_list$SE
  
  # Set random seed to allow for repeatability
  set.seed(1)
  
  # Perform SUSIE analysis using reference LD panel EUR
  fit = susie_rss(zscore, sqrt(LDmatrix), n=n)
  
  # Plot posterior inclusion probability
  Plot_PIP = susie_plot(fit, y='PIP', add_bar = T,add_legend = T)
  
  # identify location of causal variants
  sets = fit$sets
  
  # Plot
  #* highlight based on credible sets
  highlight = SNP_list %>%
    slice(c(unlist(fit$set$cs)))
  # plot
  plot_causal = ggplot(Locus, aes(x=BP,y=-log10(P))) +
    geom_point() +
    geom_label(
      data = highlight,
      nudge_x = 50000,
      aes(label = SNP)) +
    geom_point(data=highlight, 
               aes(BP,-log10(P)), 
               color='red',
               size=3) +
    theme_classic() +
    geom_hline(yintercept=5, linetype="dashed", 
               color = "red")
  
  return(list(fit, plot_causal, highlight, sets))
}


#########################
# 4 Perform SuSiE model #
#########################
#* Perform finemapping on GWAS studies using 

SusieLocus = get_cc(Locus, "EUR")
SusieLocus[[2]]
susie_plot(SusieLocus[[1]], y='PIP', add_bar = T, add_legend = T)

#####################################
# 5. Develop function to LD heatmap #
#####################################
#* Using SNP_list of
LDheat = function(SNP_list, LDLinkPopulation){
  
  LD = LDlinkR::LDmatrix(snps = SNP_list$SNP, pop = LDLinkPopulation, r2d = "r2",
                         token = token, file = FALSE)
  
  # filter SNP list from inputted GWAS file and remove any missing variants 
  SNP_list = SNP_list %>%
    filter(SNP %in% as.vector(LD[1]$RS_number))
  
  # Convert first row as rownames
  rownames(LD) = LD[,1]
  
  # Convert LD into a matrix and convert NAs to 0 
  LDmatrix = data.matrix(LD[-1] %>%
                           mutate(across(everything(), ~ replace_na(.x, 0))))
  
  return(LD.plot(LDmatrix, snp.positions = SNP_list$BP))
}


###########################
# 6. Produce LD heat maps #
###########################
#* this can take a while to render
LDheat(Locus, "EUR", 'Locusheatmap')


#####################
# 7. Summarise data #
#####################

# Manhatten plot showing idenitified putative causal variants
SusieLocus[[2]]

# Table giving details of putative causal variants
SusieLocus[[3]]
write_csv(SusieLocus[[3]], file = 'Results/LocusPutativeVariants.csv')

# Plot posterior probability of credible sets
susie_plot(SusieLocus[[1]], y='PIP', add_bar = T,add_legend = T)


