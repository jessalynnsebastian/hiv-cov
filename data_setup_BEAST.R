library(lubridate)
library(ape)
library(TransPhylo)
### SET DATA UP FOR BEAST ###
# In this file: reading the metadata, reading the fasta and splitting into lineages, getting dates formatted for BEAUTI2

# read metadata
meta <- read.csv("./combined_meta.csv") 
# edit fasta file to remove ones without metadata and get tip dates for BEAST
seq_2023.08 <- read.FASTA("./combined_seq.fasta")
# reorder meta to match the order of fasta file
meta <- meta[match(names(seq_2023.08), meta$seqName),]
# get indices of sequences for which we don't have a collection date or lineage or both, & get rid of them if they exist
miss_dat_or_lin <- which(is.na(meta$collection_date) | is.na(meta$Nextclade_pango))
seq_2023.08 <- seq_2023.08[-miss_dat_or_lin]
meta <- meta[-miss_dat_or_lin,]
# get date last sample for future conversion ape->TP ptree
dateLastSample <- decimal_date( max( as.Date(meta$collection_date), na.rm = T) ) # get latest date in decimal form, years
dateLastSample <- (dateLastSample - 2022)*365 # convert to days since start of 2022
# function to convert to TP tree (helper for future analyses since this file is sourced)
toTransPhylo <- function(tree){
  tree$edge.length <- 365*tree$edge.length # convert years to days
  tree <- ptreeFromPhylo(tree, dateLastSample) # convert ape phylo tree to transphylo ptree
  return(tree)
}

## Save full BITRI sequences + tip dates
if (!file.exists("combined_seq_nomiss_fasta") & !file.exists("bitri_dates_beauti.txt")){
  write.FASTA(seq_2023.08, "combined_seq_nomiss.fasta")
  dates_for_beauti <- meta[,c("seqName", "collection_date")]
  for ( i in 1:nrow(dates_for_beauti) ){
    cat(dates_for_beauti[i,1], "\t", dates_for_beauti[i,2], "\n", sep = "", file = paste("bitri_dates_beauti.txt", sep=""), append = T)
  }
}

## Next splitting by lineage. 
# create new column of simplified pango lineage
meta$pango_simplified <- sub("\\.(.*)","",meta$Nextclade_pango)

# Split the fasta file into 5. Already ordered meta and fasta to be in the same order so can just use indices directly
BE_fasta <- seq_2023.08[meta$pango_simplified == "BE"]
XBB_fasta <- seq_2023.08[meta$pango_simplified == "XBB"]
BA_fasta <- seq_2023.08[meta$pango_simplified == "BA"]
BQ_fasta <- seq_2023.08[meta$pango_simplified == "BQ"]
FN_fasta <- seq_2023.08[meta$pango_simplified == "FN"]

# save fasta files (if not already saved)
if (!file.exists("./BE_beast2/BE_sequences.fasta")) write.FASTA(BE_fasta, "./BE_beast2/BE_sequences.fasta")
if (!file.exists("./XBB_beast2/XBB_sequences.fasta")) write.FASTA(XBB_fasta, "./XBB_beast2/XBB_sequences.fasta")
if (!file.exists("./BA_beast2/BA_sequences.fasta")) write.FASTA(BA_fasta, "./BA_beast2/BA_sequences.fasta")
if (!file.exists("./BQ_beast2/BQ_sequences.fasta")) write.FASTA(BQ_fasta, "./BQ_beast2/BQ_sequences.fasta")
if (!file.exists("./FN_beast2/FN_sequences.fasta")) write.FASTA(FN_fasta, "./FN_beast2/FN_sequences.fasta")

## Then create tip dates in format for BEAST for each  (if not already created)
idxs <- list(BE_idx = meta$pango_simplified == "BE",
             XBB_idx = meta$pango_simplified == "XBB",
             BA_idx = meta$pango_simplified == "BA",
             BQ_idx = meta$pango_simplified == "BQ",
             FN_idx = meta$pango_simplified == "FN")
names <- c("BE", "XBB", "BA","BQ","FN")
for ( j in 1:5 ) {
  path <- paste("./", names[j], "_beast2/", names[j], "_dates_beauti.txt", sep="")
  if (!file.exists(path)){}
    meta_lineage <- meta[idxs[[j]],]
    dates_for_beauti <- meta_lineage[,c("seqName", "collection_date")]
    for ( i in 1:nrow(dates_for_beauti) ){
      cat(dates_for_beauti[i,1], "\t", dates_for_beauti[i,2], "\n", sep = "", file = path, append = T)
    }
  }
}

