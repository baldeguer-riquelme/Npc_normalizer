#!/usr/bin/env Rscript

# Author: Borja Aldeguer-Riquelme
# Email: briquelme3@gatech.edu

# Module 1. Calculate coverage for each sample
module_1 <- function(){
  npo_files<-list.files(pattern="*.npo")
  
  if (length(npo_files) == 0){
    stop(paste("\n[", Sys.time(), "] ERROR! No .npo files found in the working directory. Please, put your .npo files in the working directory", "\n"))
  }
  
  nonpareil_coverage <- data.frame()
  for (i in npo_files){
    curva<-Nonpareil.curve(i)
    label = curva@label
    C = curva@C
    number_reads <- curva@R
    nonpareil_coverage <- rbind(nonpareil_coverage, c(label, C, number_reads))
  }
  
  colnames(nonpareil_coverage) <- c("Metagenome","Coverage","Number_reads")
  nonpareil_coverage$Coverage <- as.numeric(nonpareil_coverage$Coverage)
  
  return(nonpareil_coverage)
}

# Module 2. Parse input files according to data format
parse_TAD_GEQ <- function(abund_table_3col, feature_metadata, metag_metadata){
  # Merge abundance table and feature table
  var = colnames(feature_metadata)[1]
  abund_table_3col_tad_geq <- inner_join(abund_table_3col, feature_metadata, by=c(Feature = var))
  
  # Merge abundance and feature metadata with metagenome metadata 
  abund_table_3col_tad_geq <- inner_join(abund_table_3col_tad_geq, metag_metadata, by="Metagenome")
  
  # Calculate sequencing depth
  abund_table_3col_tad_geq$SeqDepth <- abund_table_3col_tad_geq$Abund * abund_table_3col_tad_geq$GEQ
  
  return(abund_table_3col_tad_geq)
}


parse_RPKM <- function(abund_table_3col, feature_metadata, metag_metadata, nonpareil_coverage){
  # Merge abundance table and feature table
  var = colnames(feature_metadata)[1]
  abund_table_3col_RPKM <- inner_join(abund_table_3col, feature_metadata, by=c(Feature = var))
  
  # Merge abundance and feature metadata with metagenome metadata 
  #abund_table_3col_RPKM <- inner_join(abund_table_3col_RPKM, metag_metadata, by="Metagenome")
  
  #Merge abundance and metadata with nonpareil coverage data
  abund_table_3col_RPKM <- inner_join(abund_table_3col_RPKM, nonpareil_coverage, by="Metagenome")
  abund_table_3col_RPKM$Number_reads <- as.numeric(abund_table_3col_RPKM$Number_reads)
                                      
  # Calculate sequencing depth
  abund_table_3col_RPKM$Recruited_Reads <- abund_table_3col_RPKM$Abund * (abund_table_3col_RPKM$Number_reads / 1e6) * (abund_table_3col_RPKM$Length / 1e3)
  
  return(abund_table_3col_RPKM)
}


parse_reads <- function(abund_table_3col, feature_metadata, metag_metadata){
  colnames(abund_table_3col) <- c("Feature","Metagenome","Reads")
  
  # Merge abundance table and feature table
  var = colnames(feature_metadata)[1]
  abund_table_3col_reads <- inner_join(abund_table_3col, feature_metadata, by=c(Feature = var))
  
  # Merge abundance and feature metadata with metagenome metadata 
  #abund_table_3col_reads <- inner_join(abund_table_3col_reads, metag_metadata, by="Metagenome")
  
  return(abund_table_3col_reads)
}

# Integrate all functions to parse input table
module_2 <- function(abund_table_3col, feature_metadata, metag_metadata, nonpareil_coverage, input_format){
  if (input_format == "TAD_GEQ"){
    abund_table_3col_parsed <- parse_TAD_GEQ(abund_table_3col, feature_metadata, metag_metadata)
  } else if (input_format == "RPKM"){
    abund_table_3col_parsed <- parse_RPKM(abund_table_3col, feature_metadata, metag_metadata, nonpareil_coverage)
  } else if (input_format == "reads"){
    abund_table_3col_parsed <- parse_reads(abund_table_3col, feature_metadata, metag_metadata)
  }
  
  return(abund_table_3col_parsed)
}

# Module 3. Estimate recruited reads
#### Function to integrate all estimation steps
estimate <- function(sample, min_avg_coverage, abund_table_3col, input_format){
  # Get data from .npo 
  input <- paste(sample, ".npo", sep='')
  target_coverage <- min_avg_coverage
  df_reads_coverage <- data.frame()
  
  curva<-Nonpareil.curve(input)
  
  read_length <- curva@L
  number_reads <- curva@R
  
  # Calculate sequencing effort needed to reach the min_avg_coverage
  seq_effort = summary.Nonpareil.Curve(curva)["LR"]
  coverage = summary.Nonpareil.Curve(curva)["C"]
  while (coverage > target_coverage){
    seq_effort = seq_effort - 5e5 # Decrease 0.5Mb per cycle
    coverage = predict.Nonpareil.Curve(curva, lr=seq_effort)
  }
  seq_effort_round = round((seq_effort / read_length), 0)
  
  # Calculate the fraction it represent from the total
  Fraction <- seq_effort_round / number_reads 
  
  # Subsample table to get only rows with the sample of interest
  sub_samp <- abund_table_3col %>% filter(Metagenome == sample)
  
  # Estimate recruited reads at target coverage
  if (input_format == "TAD_GEQ"){
    sub_samp$SeqDepth_estim <- sub_samp$SeqDepth * Fraction
    sub_samp$GEQ_estim <- sub_samp$GEQ * Fraction
    
    sub_samp$SeqDepth_estim[sub_samp$SeqDepth_estim <= 0.1] <- 0 # Replace values below 0.1 TAD by 0
    
    sub_samp$RelAbund_estim <- sub_samp$SeqDepth_estim / sub_samp$GEQ_estim
    
  } else if (input_format == "RPKM"){
    
    sub_samp$Recruit_reads_estim <- sub_samp$Recruited_Reads * Fraction
    
    sub_samp$Recruit_reads_estim[((sub_samp$Recruit_reads_estim * read_length)/sub_samp$Length) <= 0.1] <- 0 # Replace recruited reads that cover less than 10% of the target for 0
    
    sub_samp$RelAbund_estim <- sub_samp$Recruit_reads_estim / (sub_samp$Length / 1e3) / (seq_effort_round / 1e6)
    
  }  else if (input_format == "reads"){
    
    sub_samp$RelAbund_estim <- round(sub_samp$Reads * Fraction, 0)
    
  }
  
  return(sub_samp)
}


#### Function to integrate module 3
module_3 <- function(abund_table_3col, min_avg_coverage, input_format){
  # Iterate through each sample subject to normalization
  df_estimated <- data.frame()
  for (samp in unique(abund_table_3col$Metagenome)){
    res <- estimate(samp, min_avg_coverage, abund_table_3col, input_format)
    df_estimated <- rbind(df_estimated, res)
  }
  
  # Save matrix
  Norm_matrix <- acast(formula = Feature ~ Metagenome, value.var = "RelAbund_estim", data=df_estimated) %>% as.data.frame()
  out_name <- paste("Normalized_abundance_Npc_", round(min_avg_coverage,2), ".txt",sep="")
  cat(paste("\n[ ", Sys.time(), " ] Saving normalized matrix to file: ", out_name, "\n", sep=""))
  write.table(Norm_matrix, file=out_name, sep="\t", col.names = T, row.names=T, dec=".", quote=F)
  
  return(df_estimated)
}


# Main function
library(argparser)

# Create a parser
p <- arg_parser("Estimate relative abundance (TAD/GEQ, RPKM or number of recruited reads) at a given nonpareil coverage")

# Add command line arguments
p <- add_argument(p, "--abund_matrix", help = "Relative abundance table in a matrix shape table with samples in columns and features (MAGs, genes) in rows. Relative abundance can be number of reads, RPKM or TAD/GEQ. Default is TAD/GEQ, specify a different format with argument --in_format. Required", type = "character")
p <- add_argument(p, "--in_format", help = "Input data format. Can be one of: 'reads', 'RPKM', 'TAD_GEQ'. Optional.", default='TAD_GEQ')
p <- add_argument(p, "--feature_md", help = "MAG/gene metadata. At least 2 columns are needed and should be named as 'MAG', 'Length'. Required", type = "character")
p <- add_argument(p, "--metagenome_md", help = "Metagenome metadata. At least 2 columns are needed and should be named as 'Metagenome', 'GEQ'. Only required for TAD/GEQ data.", type = "character")
p <- add_argument(p, "--nonp_cov", help = "Nonpareil coverage value to estimate relative abundance. [ 0 - 1 ]. Optional.", type = "numeric")


# Parse the command line arguments
argv <- parse_args(p)

input_abund <- argv$abund_matrix # Abundance (TAD80/GEQ) matrix table with samples in columns and features (MAGs, genes) in rows
input_format <- argv$in_format 
input_f_md <- argv$feature_md # MAG/gene metadata. At least 3 columns are needed and should be named as "MAG", "Length", "Any factor we are interested in". Any additional factors of interest can be added in next columns.
input_mg_md <- argv$metagenome_md # Metagenome metadata. At least 3 columns are needed and should be named as "Metagenome", "GEQ", "Any factor we are interested in". Any additional factors of interest can be added in next columns.
cov_to_norm <- argv$nonp_cov

# Load libraries
library(tidyverse)
library(Nonpareil)
library(reshape2)

# 0. Read tables
cat(paste("\n[", Sys.time(), "] Reading files...", "\n"))

abund_table <- read.table(input_abund, header=T, dec=".", sep="\t")
abund_table_3col <- setNames(melt(abund_table), c("Feature","Metagenome","Abund"))

feature_metadata <- read.table(input_f_md, header=T, dec=".", sep="\t")

if (is.na(input_mg_md) == FALSE){
  metag_metadata <- read.table(input_mg_md, header=T, dec=".", sep="\t",)
} 


# 1. Determine factors (Only for the future implementations, not useful now)
#list_feature_factors <- colnames(feature_metadata)[3:ncol(feature_metadata)]
#list_metag_factors <- colnames(metag_metadata)[3:ncol(metag_metadata)]


# 2. call module 1 (calculates the nonpareil coverage for each metagenome using .npo files)
cat(paste("\n[ ", Sys.time(), " ] Calculating nonpareil coverage from .npo files...", "\n", sep=""))
nonpareil_coverage <- module_1()


# 3. Check if nonp_cov was provided by user. If so, use this value for normalization. Otherwise, use the minimum nonpareil coverage
if (cov_to_norm != ""){
  if (class(cov_to_norm) == "numeric"){
    min_avg_coverage <- cov_to_norm
    cat(paste("\n[ ", Sys.time(), " ] Using nonpareil coverage value provided by user (", min_avg_coverage, ") to normalize relative abundance...", "\n", sep=""))
    if (input_format == "RPKM" & min_avg_coverage <= 0.3){
      cat(paste("\n[ ", Sys.time(), " ] WARNING! RPKM values from metagenomes with a nonpareil coverage =< 0.3 are unreliable. Results will likely be biased", sep=""))
    }
  } else {
    stop(paste("\n[ ", Sys.time()," ] ERROR! ", cov_to_norm, " is not numeric. Please, provide a value between 0 to 1", sep = ""))
  }
} else {
  min_avg_coverage <- min(nonpareil_coverage$Coverage)
  cat(paste("\n[ ", Sys.time(), " ] Using minimum nonpareil coverage value (", round(min_avg_coverage,2), ") to normalize relative abundance...", "\n", sep=""))
  if (input_format == "RPKM" & min_avg_coverage <= 0.3){
    cat(paste("\n[ ", Sys.time(), " ] WARNING! RPKM values from metagenomes with a nonpareil coverage =< 0.3 are unreliable. Results will likely be biased", sep=""))
  }
}


# 4. Parse data to get number of reads recruited for each type of data (TAD/GEQ, RPKM or reads)
abund_table_3col_parsed <- module_2(abund_table_3col, feature_metadata, metag_metadata, nonpareil_coverage, input_format)


# 5. Estimated relative abundance at given nonpareil coverage and save table
cat(paste("\n[ ", Sys.time(), " ] Normalizing data...", "\n", sep=""))
df_estimated <- module_3(abund_table_3col_parsed, min_avg_coverage, input_format)

