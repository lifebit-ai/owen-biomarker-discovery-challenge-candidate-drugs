# Methodology and Protocol Design for Identifying Candidate Drugs for Asthma

This repository contains the complete methodology and scripts for identifying potential drug targets for Asthma. 
The analysis pipeline is based on data from Genome-Wide Association Studies (GWAS), eQTL data, linkage disequilibrium data, 
and gene set enrichment analyses. The primary programming language used for this analysis is R.

# Requirements
This repository does not include the asthma summary statistics file required to replicate the analysis. 
This can be manually done by downloading data from the GWAS Catalog online database (refer to: https://www.ebi.ac.uk/gwas/studies/GCST006862). This needs to be saved to the 'raw data' folder 

# Section A:  Methodology Protocol Design
For protocol Design open word document:  

### 'Methodology and protocol design for identifying candidate drugs for Asthma.docx'

# Section B: Solution Implementation
All coding scripts used in solution implementation are detailed in Section A.  
This includes the following R scripts:

### 1. GWAStoGRange.r:

The imported summary statistics data file is standardised into a GenomicRange format.

### 2. DefineLocus.r

This visualised summary statistics data and defines locus on chromosome 17

### 3. statisticalFineMap.r

This performs a statistical fine-mapping on identified locus on chromosome 17, although this should be replicated for all loci.

### 4. eQTL.R

In the example code the following eQTL profiles were examined: IKZF3, GSDMB, ORLD3, GSDMA, PGAP3, MSL1, PNMT. These were explored in
 different single-cell tissues including: CD4_T-cell, Macrophage, lung, esophagus_mucosa, esophagus_muscularis. By examining the effects
 of these putative causal variants on gene expression in various tissues, we can gain further insights into the genetic underpinnings of Asthma. 

### 5. coloc.r

In the provided example code, colocalization was examined between the previously identified locus on chromosome 17 and eQTL data for GSDMB in lung tissue.
While the script is incomplete due to the inability to procure an LD matrix for eQTL data, the code structure remains insightful for understanding the overall process.
This colocalization analysis should be repeated for each eQTL profile to validate shared causal variant between putative causal asthma variant and gene expression.

## Visualising solution
You can view selected sections of the following files: GWAStoGRange.r, DefineLocus.r, and StatisticalFineMap.r in the document titled SolutionImplementation.html.

