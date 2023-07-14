###############
# Description #
###############
#* Creating manahatten plots from summary statistic files which have been
#* previously been tidied in GWAStoGrange.R. This allows user to zoom into a 
#* region of interest to allow infer locus.
#* 
#* This involves 
#* 1. importing genomic ranges 
#* 2. Create Manhatten plots
#* 3. Zoom into region of interest
#* 4. Export locus

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
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(karyoploteR)
require(sqldf)

##################################
# 1.Read in Gr_ranges if created #
##################################
#* Gr_ranges of chrm6 should have already been created using GWAStoGenomicRangesChrm6.R.
#* then these can be imported back in rather than needing to import large files
#* in above script.

# 1.1 Create function to convert datatable to genomic range ------------------

gRange = function(x){
  grange = makeGRangesFromDataFrame(x,
                                    seqnames.field="seqnames",
                                    start.field="start",
                                    end.field="end")
  grange$pval = x$P
  grange$SE = x$SE
  grange$BETA = x$BETA
  names(grange) = x$SNP
  return(grange)
}

# 1.2 Run function  -----------------------------------------------

GWAS = read_tsv('Grange/Granges_GWASGRCh37.tsv')
GWAS = gRange(GWAS)

#######################################
# 2. Visualise data as Manhatten plot #
#######################################
#* Create Manhatten plots using the GWAS information above. This process involves
#* i) Create colour range for GWAS data points, ii) plotting data GWAS as Manhatten.

# 2.1. Add colour scheme to data points ----------------------------------------
#* Create function
points.col = function(grRange, startCol, endCol){
  pointsCol = -log10(grRange$pval)
  pointsCol = colByValue(pointsCol, colors=c(startCol, endCol))
  return(pointsCol)
}

#* GWAS data points with grey/red colour scheme
GWAS.col = points.col(GWAS, "#BBBBBB00", "red")

# 2.2. Create reference plot 1-----------------------------------------------
kp <- plotKaryotype(plot.type=4, genome = 'hg19')
kpAddBaseNumbers(kp)

# 2.3. Create Manhatten plot-----------------------------------------------------------

kpAxis(kp, ymin=0, ymax= max(-log10(GWAS$pval)), cex = 0.5)
kpAddLabels(kp, labels = "-log10(P)", srt=90, pos=3, cex = .45)
kp = kpPlotManhattan(kp, GWAS, points.col= GWAS.col)


###################################
# 3. Zoom into region of interest #
###################################
#* From the previous manhatten plot we saw the largest effect size on chr17. We
#* can explore this further by zooming into this region. 

# Create function to find top SNP -----------------------------------------
#* requires grange file and chromosome as inputs
topSNP = function(Gr_range, CHR){
  # Filter by Chromosome 
  Gr_range = Gr_range[seqnames(Gr_range) == CHR]
  # identify top SNP in specified chromosome
  SNP = Gr_range[Gr_range@ranges@start>min(start(Gr_range)) & Gr_range@ranges@start<max(start(Gr_range))]
  topSNP = SNP[which.min(SNP$pval)]
  # Get y value (-log10(pval))
  topSNP$y = -log10(topSNP$pval)
  # remove top SNPs outside -log10 p-value threshold
  if(!is.null(which(topSNP$y< -log10(5e-8))) == TRUE){
    topSNP[which(topSNP$y< -log10(5e-8))] = NULL
  }
  return(topSNP)
}
# Identify peak SNP
topSNP(GWAS, 'chr17')

# Zoom into peak SNP region on Chrm17
#* Peak SNP was found at 39905964 so we will zoom into +/- ~1000000
kp <- plotKaryotype(plot.type=4, genome = 'hg19', zoom = 'chr17:39000000-41000000')
kpAddBaseNumbers(kp, tick.dist = 1e6/4)

# Create Manhatten plot-----------------------------------------------------------

kpAxis(kp, ymin=0, ymax= max(-log10(GWAS$pval)), cex = 0.5, r0 = 0.25)
kpAddLabels(kp, labels = "-log10(P)", srt=90, pos=3, cex = .45, r0 = 0.25)
kp = kpPlotManhattan(kp, GWAS, points.col= GWAS.col, r0 = 0.25)

# Add curated gene track --------------------------------------------
genes.data <- makeGenesDataFromTxDb(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot = kp)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpPlotGenes(kp, data=genes.data, add.transcript.names = FALSE, r1=0.25, cex=0.4, gene.name.position = "left", gene.name.cex = 0.4)
kpAddLabels(kp, labels = "Gene Track", srt=90, pos=3, r0=autotrack(1,4), cex = 0.7)

# Define locus
Locus = GRanges(seqnames = 'chr17',
                     ranges = IRanges(start = 39500000, end = 40250000))
#* Add to track
kpRect(kp, data=Locus, y0=0, y1=1, col=NA, border="red", lwd=3)


############################
# Export locus of interest #
############################

# Create an IRanges object with range of interest
range_of_interest = GRanges(
  seqnames = "chr17",
  ranges = IRanges(start = 39500000, end = 40250000)
)

# Filter the original GWAS to include only specified Locus
Locus_Gr <- subsetByOverlaps(GWAS, range_of_interest)
# specifiy SNP
Locus_Gr$SNP = names(Locus_Gr)

# Export
write.table(as.data.frame(Locus_Gr),
            file="Grange/Locus.tsv", 
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", na="")

