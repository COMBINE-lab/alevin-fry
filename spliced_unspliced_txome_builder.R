#!/usr/bin/env Rscript

#-----------------------------------------------------------------#
# Load packages
#-----------------------------------------------------------------#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("eisaR", quietly = TRUE)) {
  BiocManager::install("eisaR")
}
if(!requireNamespace("tximeta", quietly = TRUE)) {
  BiocManager::install("tximeta")
}
suppressPackageStartupMessages({
  library(eisaR)
  library(BSgenome)
})
#-----------------------------------------------------------------#
# Initialization
#-----------------------------------------------------------------#

help_info = "
This Rscript is used for building expanded transcriptome containing the spliced and unspliced version of genes.
-----------------------------------------------------------
Inputs:
  1. --gtf: a gtf/gtf(.gz) file
  2. --fa: a genome.fa(.gz) file
  3. (optional) --outdir: an output directory

output:
An expanded transcriptome contains:
  1. annotation.expanded.tx2gene.tsv: Transcript names with corresponding gene names
  2. annotatoion.expanded.gtf: Each gene has an unspliced version, which is the genome range from the 5' end of the first exon to the 3' end of the last exon
  3. annotation.expanded.fa: The corresponding sequences of the terms in annotatoion.expanded.gtf
Note: 
  1. The expanded transcriptome is named using the name of the input gtf file
  2. The correct input arg example: --fa=genome.fa
  3. If the outdir is not specified, a new directory will be created in the working directory 
  4. If some strange errors show up, please check the version of the packages, this script requires the same environment as `eisaR` R package.
  5. This script is based on `eisaR` and `tximeta` R packages.
  
"

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  # If there is only one arg, it should be --help
  fa = NaN
  gtf = NaN
  help_mode = FALSE
  outdir = file.path(getwd(),"expanded.txome")
  if (length(args)==1 & args[1] =="--help") {
    message(help_info)
    help_mode = TRUE
  } else {
    for (arg in args) {
      if (startsWith(arg,"--fa")) {
        fa = strsplit(arg, "--fa=")[[1]][2]
      } else if (startsWith(arg,"--gtf")) {
        gtf = strsplit(arg, "--gtf=")[[1]][2]
      }  else if (startsWith(arg,"--outdir")) {
        outdir = strsplit(arg, "--outdir=")[[1]][2]
        if(!endsWith(outdir,"/")|!endsWith(outdir,"\\")){
          outdir = file.path(outdir,"")
        }
      }
    }
  }
  if ((is.na(fa)|is.na(gtf))&!help_mode){
    stop("Please use --help to check expected inputs.", call.=FALSE)
  } else if(!help_mode) {
    process(gtf=gtf, fa=fa, outdir=outdir)
  }
}

#-----------------------------------------------------------------#
# main step
#-----------------------------------------------------------------#

process <- function(gtf, fa, outdir) {
  # setup output file name
  dir.create(file.path(outdir), showWarnings = FALSE)
  out_fa=file.path(outdir,"annotation.expanded.fa")
  out_gtf=file.path(outdir,"annotation.expanded.gtf")
  out_tx2gene=file.path(outdir,"annotation.expanded.tx2gene.tsv")
  
  # Initialization
  featureType = c("spliced", "unspliced")
  suffixes <- c(unspliced = "-U")
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
  grlfull <- GenomicRanges::GRangesList()
  featurelist <- list()
  ebt <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  t2g <- AnnotationDbi::select(txdb, keys = names(ebt), keytype = "TXNAME", 
                               columns = "GENEID")
  e2 <- BiocGenerics::unlist(ebt)
  e2$transcript_id <- names(e2)
  e2$gene_id = t2g$GENEID[match(e2$transcript_id, t2g$TXNAME)]
  e2$exon_id <- e2$exon_name
  e2$exon_name <- NULL
  e2$type <- "exon"
  names(e2) <- NULL
  S4Vectors::mcols(e2) <- S4Vectors::mcols(e2)[, c("exon_id", "exon_rank", 
                                        "transcript_id", "gene_id", "type")]
  ebt <- BiocGenerics::relist(e2, ebt)
  ebg <- GenomicFeatures::exonsBy(txdb, by = "gene")
  corrtx <- data.frame(spliced = unique(names(ebt)), stringsAsFactors = FALSE)
  corrtx$unspliced <- paste0(corrtx$spliced, suffixes["unspliced"])
  corrgene <- data.frame(spliced = unique(unlist(ebt)$gene_id),
                         stringsAsFactors = FALSE)
  corrgene$unspliced <- paste0(corrgene$spliced, suffixes["unspliced"])
  
  # process spliced txps
  message("Extracting spliced transcript features")
  featurelist$spliced <- names(ebt)
  grlfull <- c(grlfull, ebt)
  
  #process unspliced txps
  message("Extracting unspliced transcript features")
  ebtr <- range(ebg)
  e2 <- BiocGenerics::unlist(ebtr)
  
  e2$exon_rank <- 1L
  e2$transcript_id <- names(e2)
  e2$gene_id = names(e2)
  
  e2$type <- "exon"
  e2$transcript_id <- paste0(e2$transcript_id, suffixes["unspliced"])
  e2$exon_id <- e2$transcript_id
  names(e2) <- NULL
  S4Vectors::mcols(e2) <- S4Vectors::mcols(e2)[, c("exon_id", "exon_rank", 
                                        "transcript_id", "gene_id", "type")]
  ebtr <- BiocGenerics::relist(e2, ebtr)
  names(ebtr) <- paste0(names(ebtr), suffixes["unspliced"])
  featurelist$unspliced <- names(ebtr)
  grlfull <- c(grlfull, ebtr)
  
  # post processing
  S4Vectors::metadata(grlfull)$corrtx <- corrtx[, colnames(corrtx) %in%
                                                  featureType, drop = FALSE]
  S4Vectors::metadata(grlfull)$corrgene <- corrgene[, colnames(corrgene) %in%
                                                      featureType, drop = FALSE]
  S4Vectors::metadata(grlfull)$featurelist <- featurelist
  grl = grlfull
  
  # rm(corrgene,corrtx,e2,ebg,ebt,ebtr,featurelist,grlfull,spliced_txp_length,t2g,unspliced_txp_length)
  
  # writing outputs
  message("Reading genome...")
  
  genome <- Biostrings::readDNAStringSet(
    fa
  )
  names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
  
  message("Extracting transcript sequences...")
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = genome, 
    transcripts = grl
  )
  
  message("Writing transcript sequences...")
  Biostrings::writeXStringSet(
    seqs, out_fa
  )

  message("Writing gtf file...")
  eisaR::exportToGtf(
    grl, 
    filepath = out_gtf
  )
  
  message("Writing tx2gene file...")
  df <- eisaR::getTx2Gene(
    grl, filepath = out_tx2gene
    )
  
}

main()