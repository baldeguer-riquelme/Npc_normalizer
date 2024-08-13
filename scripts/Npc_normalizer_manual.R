#!/usr/bin/env Rscript

# Author: Borja Aldeguer-Riquelme
# Email: briquelme3@gatech.edu


# Main function
library(argparser)
library(Nonpareil)


# Create a parser
p <- arg_parser("Subsample a metagenome to a given nonpareil coverage")

# Add command line arguments
p <- add_argument(p, "--m", help = "Path to metagenome reads in fastq or fasta format. Allowed extensions: fastq, fastq.gz, fastq.bz2, fasta, fasta.gz, fasta.bz2. Required.", type = "character")
p <- add_argument(p, "--cov", help = "Target Nonpareil coverage value. From 0 to 1. Required.", type = "double")
p <- add_argument(p, "--m2", help = "For pair end reads, path to reverse reads in fastq or fasta format. Allowed extensions: fastq, fastq.gz, fastq.bz2, fasta, fasta.gz, fasta.bz2. Optional.", type = "character")
p <- add_argument(p, "--npo", help = "Path to .npo files. By default, will look for files with the same name as in --m but with extension .npo. Optional unless npo files have a different name structure.", type = "character")

# Parse the command line arguments
argv <- parse_args(p)



# 0. Check --m and --cov arguments exists
if (is.na(argv$m) == TRUE){
  cat("ERROR. Argument --m is required. Please provide the path to the metagenome reads to subsample.")
  stop()
}

if (is.na(argv$cov) == TRUE){
  cat("ERROR. Argument --cov is required. Please provide the Nonpareil coverage to subsample your metagenomes.")
  stop()
}




# 1. Find npo file using base name
base_name <- sub("\\.fast.*$", "", argv$m)
extension <- sub(".*\\.fast", ".fast", argv$m)

if (is.na(argv$npo) == TRUE){
  npo <- paste(base_name, ".npo", sep="")
} else {
  npo <- argv$npo
}




# 2. Calculate number of reads in subsampled metagenomes
curva <- Nonpareil.curve(npo)
read_total <- curva@R
read_length <- curva@L

cat(paste("\n[", Sys.time(), "] Calculating number of reads for Nonpareil coverage ", argv$cov, "...", "\n", sep=""))

seq_effort = summary.Nonpareil.Curve(curva)["LR"]
coverage = summary.Nonpareil.Curve(curva)["C"]
decrease_ratio = seq_effort*0.0001
while (coverage > argv$cov){
  seq_effort = seq_effort - decrease_ratio # Decrease 0.01% reads per cycle
  coverage = predict.Nonpareil.Curve(curva, lr=seq_effort)
}

read_subsampling = round((seq_effort / read_length), 0)

cat(paste("\n\nInitial number of reads: ", read_total, "\nNumber of reads for Nonpareil coverage ", argv$cov, ": ", read_subsampling,"\n\n", sep=""))




# 3. Subsample metagenomes
if(is.na(argv$m2) == TRUE){
  
  cat(paste("\n[", Sys.time(), "] Subsampling metagenome ", argv$m, "...", "\n\n", sep=""))
  
  command <- paste("reformat.sh in=", argv$m, 
                   " out=", base_name, "_Npc_", argv$cov, extension, 
                   " samplereadstarget=", read_subsampling,
                   sep="")
  
  system(command)
  
} else {  
  
  cat(paste("\n[", Sys.time(), "] Subsampling metagenomes ", argv$m, " ", argv$m2, "...", "\n\n", sep=""))
  
  base_name1 <- sub("\\.fast.*$", "", argv$m)
  extension1 <- sub(".*\\.fast", ".fast", argv$m)
  
  base_name2 <- sub("\\.fast.*$", "", argv$m2)
  extension2 <- sub(".*\\.fast", ".fast", argv$m2)
  
  command <- paste("reformat.sh in=", argv$m,
                   " in2=", argv$m2,
                   " out=", base_name1, "_Npc_", argv$cov, extension1,
                   " out2=", base_name2, "_Npc_", argv$cov, extension2,
                   " samplereadstarget=", read_subsampling,
                   sep="")
  
  system(command)
  
}
