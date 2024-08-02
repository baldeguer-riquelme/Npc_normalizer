#!/usr/bin/env Rscript

# Author: Borja Aldeguer-Riquelme
# Email: briquelme3@gatech.edu

# Module 1. Calculate coverage for each sample
module_1 <- function(){
  npo_files<-list.files(pattern="*.npo")
  
  nonpareil_coverage <- data.frame()
  for (i in npo_files){
    curva<-Nonpareil.curve(i)
    C = curva@C
    label = curva@label
    diversity <- curva@diversity
    nonpareil_coverage <- rbind(nonpareil_coverage, c(label, C, diversity))
  }
  
  colnames(nonpareil_coverage) <- c("Metagenome","Coverage", "Diversity")
  nonpareil_coverage$Coverage <- as.numeric(nonpareil_coverage$Coverage)
  nonpareil_coverage$Diversity <- as.numeric(nonpareil_coverage$Diversity)
  
  return(nonpareil_coverage)
}


# Module 2. Estimate TAD80 values for the group of metagenomes with the highest average value
#### Function to estimate sequencing depth
estimate_SD_GEQ = function(sample_to_norm_metadata, x, df_reads_coverage, read_length){
  new_df <- sample_to_norm_metadata[,c("Feature", "Metagenome", "Length")]
  
  cov <- df_reads_coverage[x,]$Coverage
  fraction <- df_reads_coverage[x,]$Fraction
  new_df['Coverage'] <- cov
  
  new_df['Recruited_reads'] <- sample_to_norm_metadata$Recruited_reads * fraction
  estim_recruit_bases <- new_df['Recruited_reads'] * read_length
  
  new_df['SeqDepth'] <- round((estim_recruit_bases  / sample_to_norm_metadata$Length),3)
  new_df['GEQ'] <- sample_to_norm_metadata$GEQ * fraction
  new_df['RA'] <- new_df['SeqDepth'] / new_df['GEQ']
  
  # return 
  return(new_df)
}


#### Function to integrate all estimation steps
estimate <- function(sample_to_norm_metadata, median_cov, sample, npo_path){
  input <- paste(npo_path, sample, ".npo", sep='')
  df_reads_coverage <- data.frame()
  
  curva<-Nonpareil.curve(input)
  read_length <- curva@L
  reads_total <- curva@R
  
  seq_effort = summary.Nonpareil.Curve(curva)["LR"]
  coverage = summary.Nonpareil.Curve(curva)["C"]
  decrease_ratio = seq_effort*0.0001
  for (target_coverage in seq(median_cov, 0.1, by=-0.01)){
    while (coverage > target_coverage){
      seq_effort = seq_effort - decrease_ratio # Decrease 0.01% per cycle
      coverage = predict.Nonpareil.Curve(curva, lr=seq_effort)
    }
    reads_subsampling <- round((seq_effort / read_length), 0)
    df_reads_coverage <- rbind(df_reads_coverage, c(sample, target_coverage, reads_subsampling))
  }
  
  colnames(df_reads_coverage) <- c("Sample", "Coverage", "Reads_subsampling")
  df_reads_coverage$Reads_subsampling <- as.numeric(df_reads_coverage$Reads_subsampling)
  
  df_reads_coverage$Fraction <- df_reads_coverage$Reads_subsampling / reads_total
  
  # Calculate recruited reads
  sample_to_norm_metadata$Recruited_reads <- round(((sample_to_norm_metadata$Length * sample_to_norm_metadata$SeqDepth) / read_length), 0)
  sample_to_norm_metadata <- sample_to_norm_metadata[,c("Feature", "Metagenome", "Length", "Coverage", "Recruited_reads", "SeqDepth", "GEQ", "RA")]
  
  # Estimate SD
  df_estimated_SD_GEQ <- data.frame()
  for (x in 1:length(df_reads_coverage$Coverage)){
    estimated <- estimate_SD_GEQ(sample_to_norm_metadata, x, df_reads_coverage, read_length)
    df_estimated_SD_GEQ <- rbind(df_estimated_SD_GEQ, estimated)
  }
  
  # Minimum Sequencing Depth = 0.1. Consider the rest as undetected.
  df_estimated_SD_GEQ[,c("SeqDepth")][df_estimated_SD_GEQ[,c("SeqDepth")] <= 0.1] <- 0 # Replace values below 0.1 SD by 0
  df_estimated_SD_GEQ[,c("RA")][df_estimated_SD_GEQ[,c("SeqDepth")] <= 0.1] <- 0 # Replace values below 0.1 SD by 0
  
  return(df_estimated_SD_GEQ)
}

#### Function to integrate module 2
module_2 <- function(abund_table_3col_metadata, sample, median_cov, npo_path){
  # Keep only rows from the sample of interest
  sample_to_norm_metadata <- abund_table_3col_metadata %>% filter(Metagenome == sample)
  
  # Calculate SD in the original metagenome
  sample_to_norm_metadata$SeqDepth <- sample_to_norm_metadata$RA * sample_to_norm_metadata$GEQ
  
  # Iterate through each sample subject to normalization
  df_final_estimated <- estimate(sample_to_norm_metadata, median_cov, sample, npo_path)
  
  return(df_final_estimated)
}



# Module 3. Get the maximum difference in nonpareil value to get statistically significant differences
plot_avg_Npc <- function(sub_feat_level, feature_factor, feature_level, metagenome_factor, level){
  sub_feat_level_sum <- sub_feat_level %>% group_by(Metagenome, Coverage) %>% summarise(Sum=sum(RA), Num_features = sum(RA > 0))
  
  layout(mat =  matrix(c(1, 2)), heights = c(1,1)) # Heights of the two rows
  par(cex=1.5, cex.axis=.8, mar=c(2,4,2,2))
  
  boxplot(Sum ~ Coverage, data = sub_feat_level_sum, 
          main = paste(feature_factor, ": ", feature_level, "\n", metagenome_factor, ": ", level, sep = ""),
          cex.main=0.8, font.main=2,
          horizontal = F, ylab="Aggregated Relative Abundance", col = c("#0e4bac"), xaxt='n', xlab = NULL)
  
  par(cex=1.5, cex.axis=.8, mar=c(4,4,0,2))
  boxplot(Num_features ~ Coverage, data = sub_feat_level_sum, 
          las = 2,
          horizontal = F, xlab = "Nonpareil Coverage", ylab="# features", col = c("#dcc88f"))
}


stats_test <- function(feature_level, feature_factor, df_all_cov_estim_to_analyze, metagenome_factor, level){
  # Subsample df to get only the rows of the feature of interest
  sub_feat_level <- df_all_cov_estim_to_analyze %>% filter(get(feature_factor) == feature_level)
  
  # Loop over coverage values and perform t.test
  ref_coverage = max(sub_feat_level$Coverage) %>% as.double()
  num_metagenomes <- sub_feat_level %>% filter(Coverage == ref_coverage & RA > 0) %>% 
    select(Metagenome) %>% unique() %>% nrow()
  num_features_ref <- sub_feat_level %>% filter(Coverage == ref_coverage & RA > 0) %>% 
    select(Feature) %>% unique() %>% nrow()
  avg_RA <- sub_feat_level %>% filter(Coverage == ref_coverage & RA > 0) %>% summarise(Avg_RA=mean(RA)) %>% as.double()
  median_RA <- sub_feat_level %>% filter(Coverage == ref_coverage & RA > 0) %>% summarise(Avg_RA=median(RA)) %>% as.double()
  sub_feat_level$Coverage <- as.numeric(sub_feat_level$Coverage)
  
  # Plot abundance and number of features across Npc
  if (num_metagenomes > 0){
    plot_avg_Npc(sub_feat_level, feature_factor, feature_level, metagenome_factor, level)
  }
  
  # For factor levels with 2 or more metagenomes. Perform grouped by sample t.test
  if (num_metagenomes > 1){
    Diff_cov <- NA
    ratio_cov <- NA
    for (cov_level in rev(seq(0.1, ref_coverage - 0.01, by=0.01))){
      sub_feat_level_cov <- sub_feat_level %>% filter(near(Coverage, ref_coverage) | near(Coverage, cov_level)) %>%
        group_by(Metagenome, Coverage) %>% summarise(Sum=sum(RA))
      res_t.test <- t.test(Sum ~ Coverage, sub_feat_level_cov)
      if (is.nan(res_t.test$p.value) == FALSE){
        if (res_t.test$p.value < 0.05 | res_t.test$estimate[[1]] == 0){
          Diff_cov <- ref_coverage - cov_level
          ratio_cov <- Diff_cov / ref_coverage * 100
          break
        }
      }
    }
    # For factors levels with only 1 metagenome
    } else if (num_metagenomes == 1){ 
      if (num_features_ref > 1){ # If there are more than 1 features perform t.test (without grouping by sample)
        Diff_cov <- NA
        ratio_cov <- NA
        for (cov_level in rev(seq(0.1, ref_coverage - 0.01, by=0.01))){
          sub_feat_level_cov <- sub_feat_level %>% filter(near(Coverage, ref_coverage) | near(Coverage, cov_level))
          res_t.test <- t.test(RA ~ Coverage, sub_feat_level_cov)
          if (is.nan(res_t.test$p.value) == FALSE){
            if (res_t.test$p.value < 0.05 | res_t.test$estimate[[1]] == 0){
              Diff_cov <- ref_coverage - cov_level
              ratio_cov <- Diff_cov / ref_coverage * 100
              break
            }
          }
        }
      } else { # If there is only 1 feature, identify the minimum nonpareil coverage to detect that feature 
        sub_feat_level_undetected <- sub_feat_level %>% filter(RA == 0) 
        if (length(sub_feat_level_undetected$Coverage) > 0){
          cov_level <- sub_feat_level_undetected %>% summarise(Max=max(as.numeric(Coverage)) + 0.01) %>% as.double()
          Diff_cov <- ref_coverage - cov_level
          ratio_cov <- Diff_cov / ref_coverage * 100
        } else {
          Diff_cov <- NA
          ratio_cov <- NA
          cov_level <- NA
        }
      }
      # For features not present in this sample (or group of samples)
    } else if (num_metagenomes < 1){ 
      Diff_cov <- NA
      ratio_cov <- NA
      cov_level <- NA
      }
  
  # When cov_level is equal to 0.1 means that no statistically significant difference was found (0.1 is the minimum Npc tested)
  if (is.numeric(cov_level) == TRUE & near(cov_level, 0.1) == TRUE){
    cov_level <- NA 
    Diff_cov <- NA 
    ratio_cov <- NA
    } 
  
  # Return results in a list format
  return(list(feature_factor, feature_level, metagenome_factor, level, ref_coverage, cov_level, Diff_cov, ratio_cov, num_features_ref, avg_RA, median_RA))
}



# Create a parser
library(argparser)
p <- arg_parser("The script calculates ΔNpc_max")

# Add command line arguments
p <- add_argument(p, "--abund_matrix", help = "Relative abundance table in a matrix shape table with samples in columns and features (MAGs, genes) in rows. Relative abundance can be number of reads, RPKM or TAD/GEQ. Default is TAD/GEQ, specify a different format with argument --in_format. Required", type = "character")
p <- add_argument(p, "--in_format", help = "Input data format. Can be one of: 'TAD_GEQ'. Optional.", default='TAD_GEQ')
p <- add_argument(p, "--feature_md", help = "MAG/gene metadata. At least 2 columns are needed and should be named as 'MAG', 'Length'. Required.", type = "character")
p <- add_argument(p, "--metagenome_md", help = "Metagenome metadata. At least 2 columns are needed and should be named as 'Metagenome', 'GEQ'. Only required for TAD/GEQ data.", type = "character")
p <- add_argument(p, "--npo", help = "Folder containing .npo files which must have the same name as in the abundance matrix. Required.", type = "character", default="./")
p <- add_argument(p, "--out", help = "Output file name. Tab delimited.", type = "character", default="Output_Npc_max")

# Parse the command line arguments
argv <- parse_args(p)

# Load libraries
library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(tidyverse)
library(Nonpareil)
library(data.table)
library(cowplot)


# On Mac
#setwd("/Users/briquelme3/Dropbox (GaTech)/Metagenome_comparison_uneven_NonpareilCov/Metagenomic_data/tests/")


# 0. Read tables
abund_table <- read.table(argv$abund_matrix, header=T, dec=".", sep="\t")
abund_table_3col <- setNames(reshape2::melt(abund_table), c("Feature","Metagenome","RA")) %>% add_column(Coverage = "Original")
feature_metadata <- read.table(argv$feature_md, header=T, dec=".", sep="\t")
metag_metadata <- read.table(argv$metagenome_md, header=T, dec=".", sep="\t",)
npo_path <- argv$npo
out <- argv$out


# 1. call module 1 (calculates coverage for each metagenome using .npo files)
nonpareil_coverage <- module_1()


# 2. Normalize all samples using the median Npc values as the starting point
median_cov <- round(median(nonpareil_coverage$Coverage), 2)

cat(paste("[",Sys.time(),"] Normalizing abundance from ", median_cov, " to 0.1\n", sep = ""))

## Add metagenome metadata to 3col table
abund_table_3col_metadata <- inner_join(abund_table_3col, metag_metadata, by="Metagenome")

## Add feature metadata to 3col table
var = colnames(feature_metadata)[1]
abund_table_3col_metadata <- inner_join(abund_table_3col_metadata, feature_metadata, by=c("Feature" = var))

## Estimate
df_all_cov_estim <- data.frame()
for (sample in colnames(abund_table[-1])){
  cat(paste("[",Sys.time(),"] Normalizing ", sample, "\n", sep = ""))
  df_final_estimated <- module_2(abund_table_3col_metadata, sample, median_cov, npo_path)
  df_all_cov_estim <- rbind(df_all_cov_estim, df_final_estimated)
}


# 3. Calculate the maximum difference in Npc for any factor
## Determine feature factors
if (ncol(feature_metadata) > 2){
  list_feature_factors <- colnames(feature_metadata)[3:ncol(feature_metadata)]
  df_all_cov_estim <- inner_join(df_all_cov_estim, feature_metadata[,c(-2)], by=c("Feature" = var))
} else {
  list_feature_factors <- "Feature"
}

## Determine metagenome factors
if (ncol(metag_metadata) > 2){
  list_metag_factors <- colnames(metag_metadata)[3:ncol(metag_metadata)]
  df_all_cov_estim <- inner_join(df_all_cov_estim, metag_metadata[,c(-2)], by="Metagenome")
} else {
  list_metag_factors <- "Metagenome"
}

out_pdf <- paste(out, ".pdf", sep="")
pdf(out_pdf, width = 13, height = 13)
Diff_npc_values <- list()
for (feature_factor in list_feature_factors){
  for (metagenome_factor in list_metag_factors){
    levels_to_norm <- unique(metag_metadata[[metagenome_factor]])
    for (level in levels_to_norm){
      cat(paste("\n[",Sys.time(),"] Calculating ΔNpc_max for: ", feature_factor, " vs ", metagenome_factor, " ", level, "\n", sep = ""))
      
      # Get names of samples to norm
      samples_to_analyze <- metag_metadata %>% filter(get(metagenome_factor) == level) %>% select(Metagenome)
      
      # Get abundance table for samples to norm
      df_all_cov_estim_to_analyze <- df_all_cov_estim %>% filter(Metagenome %in% samples_to_analyze[,c(1)])
      
      # Get nonpareil coverage of samples to norm
      np_cov_subset <- nonpareil_coverage %>% filter(Metagenome %in% samples_to_analyze[,c(1)])
      
      # Get average nonpareil diversity
      avg_Npd <- round(mean(np_cov_subset$Diversity),2)
      
      # Get a list with the maximum difference in nonpareil values to get statistically significant differences
      list_feat_levels <- unique(df_all_cov_estim_to_analyze[[feature_factor]]) %>% as.list()
      
      list_res <- lapply(list_feat_levels, stats_test, feature_factor, df_all_cov_estim_to_analyze, metagenome_factor, level)
      for (x in seq(1,length(list_res))){
        list_res[[x]] <- append(list_res[[x]],avg_Npd)
        }
      Diff_npc_values <- append(Diff_npc_values, list_res)
      rm(list_res)
    }
  }
}
dev.off()

# Merge all lists results provided from the stats_test function to a single df
df_diff_npc_values <- rbindlist(lapply(Diff_npc_values, as.data.frame), use.names=F )
colnames(df_diff_npc_values) <- c("Feature_factor", "Feature_level", "Metagenome_factor", "Metagenome_level", "Ref_Npc", "Npc_min","ΔNpc_max (Ref_Npc-Npc_min)", "Ratio Ref_Npc/Min_Npc", "Num_features", "Avg_RA", "Median_RA", "Avg_MG_diversity")

# Remove those features that were not detected
df_diff_npc_values <- df_diff_npc_values %>% filter(Num_features > 0)

# Save results into a tab delimited table
out_tsv <- paste(out, ".tsv", sep="")
cat(paste("[",Sys.time(),"] Writing results to file: ", out_tsv, "\n\n", sep = ""))
write.table(df_diff_npc_values, out_tsv, col.names=T, row.names = F, sep = "\t", quote = F, )

