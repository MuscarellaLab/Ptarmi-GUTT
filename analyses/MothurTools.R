#                                                                              #
# MothurTools Functions Source Code                                            #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Mario Muscarella                                                 #
#                                                                              #
# Last update: 2022/08/10                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code provides functions to be used in the analysis of            #
#        16S rRNA sequence data post mothur anlaysis                           #
#                                                                              #
# Issues: Slow performance reading in OTU tables (common R issue)              #
#                                                                              #
# Recent Changes:                                                              #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. Design functions to work with shared files in memory              #
#         2. Add warnings                                                      #
#                                                                              #
################################################################################

require("reshape")||install.packages("reshape");require("reshape")

# Import OTU Site-by-Species Matrix
read.otu <- function(shared = " ", cutoff = "0.03"){
  matrix <- read.table(shared, header=T, fill=TRUE, comment.char="", sep="\t")
  matrix.cutoff <- subset(matrix, matrix$label == cutoff)
  matrix.out    <- as.matrix(matrix.cutoff[1:dim(matrix.cutoff)[1],
                                           4:(3+mean(matrix.cutoff$numOtus))])
  row.names(matrix.out) <- matrix.cutoff$Group
  return(matrix.out)
  }

# Import Taxonomy Information
read.tax <- function(taxonomy = " ", format = c("rdp", "silva"),
                     col.tax = "3", tax.levels = "6"){
  tax_raw <- read.delim(taxonomy, header = F, skip = 1)       # load genus-level data
  tax.levels <- as.numeric(tax.levels)
  col.tax <- as.numeric(col.tax)
  if (format == "rdp"){
    tax.info <- data.frame(matrix(NA, dim(tax_raw)[1], tax.levels))
    colnames(tax.info) <- c("Domain","Phylum","Class","Order","Family",
                            "Genus", "Species")[1:tax.levels]
    for (i in 1:dim(tax_raw)[1]){
      tax.split <- strsplit(as.character(tax_raw[i,col.tax]), split="\\;")[[1]]
      if (length(tax.split) < 6){
        for (j in (length(tax.split)+1):6){
             tax.split[j] = "unknown"
        }}
      tax.info[i,] <- tax.split
    }
    tax.otu.raw <- tax_raw[,1]
    tax.otu.r <- gsub(",.+$", "", tax.otu.raw)
    tax.otu <- gsub("___.+$", "_", tax.otu.r)
    #tax.otu <- sapply(strsplit(tax.otu.raw, ","), "[[", 1)
    tax <- cbind(OTU = tax.otu,tax.info)
    for (i in 2:7){
      tax[,i] <- gsub("\\(.*$", "", tax[,i])
    }
  } else {
    if (format == "silva"){
      tax.info <- data.frame(matrix(NA, dim(tax_raw)[1], tax.levels))
      colnames(tax.info) <- c("Domain","Phylum","Class","Order","Family",
                              "Genus", "Species")[1:tax.levels]
      for (i in 1:dim(tax_raw)[1]){
        tax.split <- strsplit(as.character(tax_raw[i,col.tax]), split="\\;")[[1]]
        if (length(tax.split) < 6){
          for (j in (length(tax.split)+1):6){
            tax.split[j] = "unknown"
          }}
        tax.info[i,] <- tax.split
      }
      tax.otu.raw <- tax_raw[,1]
      tax.otu.r <- gsub(",.+$", "", tax.otu.raw)
      tax.otu <- gsub("___.+$", "_", tax.otu.r)
      #tax.otu <- sapply(strsplit(tax.otu.raw, ","), "[[", 1)
      tax <- cbind(OTU = tax.otu,tax.info)
      for (i in 2:7){
        tax[,i] <- gsub("\\(.*$", "", tax[,i])
      }
      } else {
        stop("This funciton currently only works for RDP and Silva")
      }
  }
  return(tax)
}
