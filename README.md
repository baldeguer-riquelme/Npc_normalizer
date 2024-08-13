# Nonpareil coverage normalizer
Repository with scripts to normalize relative abundance (TAD/GEQ, RPKM or recruited reads) to a given Nonpareil coverage (Npc).

### Install dependencies
```
conda create -n npc_norm -c conda-forge -c bioconda r-tidyverse r-reshape2 nonpareil r-roxygen2
```
```
conda activate npc_norm
```

To install argparser:
```
# In terminal
git clone https://bitbucket.org/djhshih/argparser.git
cd argparser
R
```
```
# In R
library(roxygen2)
roxygenize()
quit()
```
```
# In terminal
R CMD INSTALL .
```

### Test installation
To verify that the script runs without errors, run a test with the files placed in the "examples" folder.

```
cp Npc_normalizer/scripts/Npc_normalizer.R Npc_normalizer/examples
cd Npc_normalizer/examples
Rscript Npc_normalizer.R -a Matrix_TAD_GEQ.txt -f Feature_metadata.txt -m Metagenome_metadata.txt -i TAD_GEQ -n 0.5 
```

The output should look like:

```
── Attaching core tidyverse packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.1.4     ✔ readr     2.1.5
✔ forcats   1.0.0     ✔ stringr   1.5.1
✔ ggplot2   3.5.1     ✔ tibble    3.2.1
✔ lubridate 1.9.3     ✔ tidyr     1.3.1
✔ purrr     1.0.2
── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

Attaching package: ‘reshape2’

The following object is masked from ‘package:tidyr’:

    smiths


[ 2024-06-18 16:22:05 ] Reading files...
Using Gene as id variables

[ 2024-06-18 16:22:05 ] Calculating nonpareil coverage from .npo files...

[ 2024-06-18 16:22:14 ] Using nonpareil coverage value provided by user (0.5) to normalize relative abundance...

[ 2024-06-18 16:22:15 ] Normalizing data...

[ 2024-06-18 16:22:47 ] Saving normalized matrix to file: Normalized_abundance_Npc_0.5.txt
```

## How to run
## 1. Check if your data requires Nonpareil coverage normalization based on ΔNpc_max
In the case that the Nonpareil coverage of your metagenome(s) is different between the factors you are interested in (e.g., treated vs control metagenomes) you should calculate ΔNpc_max to assess whether such difference could have an impact on your results. You can also calculate this metric to determine the minimum Npc needed to detect a genome/gene in one metagenome. The script Npc_max.R was developed for this purpose and you will need:

1. Abundance matrix. Should be tab delimited matrix with features (MAGs, genes) in rows and samples in columns. Only TAD/GEQ is accepted for now. We plan to incorporate RPKM in the future.
2. MAG metadata. MAG IDs and length (bp) should be provided in 2 columns named as "MAG" and "Length". Genes can be used instead of MAGs but the column header should be kept as "MAG". Feature factors (e.g., taxonomy or gene function) along with their corresponding levels (e.g., ) can be added in the third and following columns.
3. Metagenome metadata. Metagenome IDs and genome equivalents (GEQ) should be provided with at least 2 columns named as "Metagenome" and "GEQ". Metagenome factors (e.g., Antibiotic treatment) along with their corresponding levels (e.g., treated vs untreated) can be added in the third and following columns.
4. .npo files generated after running nonpareil (stand-alone) with each of your metagenomes of interest. These files should be placed on the same folder where you are running the analysis. For instructions about how to run nonpareil visit https://nonpareil.readthedocs.io/en/latest/index.html

Once you have them, you are ready to calculate ΔNpc_max running:
```
Rscript Npc_max.R -a Matrix_TAD_GEQ.txt -f Feature_metadata.txt -m Metagenome_metadata.txt --npo npo_folder/ -o output_test
```

Different statistical approaches will be used depending on the number of features and metagenomes per factor level:

|              | 1 metagenome                                        | ≥ 2 metagenomes                                         |
|--------------|-----------------------------------------------------|---------------------------------------------------------|
| 1 feature    | Minimum Npc to detect the feature in the metagenome | t.test of the aggregated relative abundances per sample |
| ≥ 2 features | t.test                                              | t.test of the aggregated relative abundances per sample |

The output will be 1) a pdf with boxplots displaying the abundance and number of features detected at each Npc and 2) a tab delimited table with the ΔNpc_max for each metagenome and feature factor. In the case that factors are not provided, the script will calculate the minimum Npc needed to detected each individual feature (genome or gene) in each metagenome. If ΔNpc_max of your feature(s) of interest is lower than the actual Npc difference between your metagenomes you will need to normalize your metagenomes to the same Npc (see section 2 below).

## 2. Normalization
In the case that your data requires Nonpareil coverage normalization, there are two options: 1) subsampling the metagenome or 2) estimate relative abundances.

### 2.1. Metagenome subsampling
To get your subsampled metagenome at a given Nonpareil coverage you will need:

1. Metagenome in FastQ/A format. Single or pair end.
2. npo files from Nonpareil analysis 

Then, you will have to choose the Nonpareil coverage value you would like to subsample your metagenome and run:
```
# For single end metagenomes
Rscript Npc_normalizer_manual.R --m metagenome.fastq --npo metagenome.npo --cov 0.6
```

```
# For pair end metagenomes
Rscript Npc_normalizer_manual.R --m metagenome_1.fastq --m2 metagenome_2.fastq --npo metagenome.npo --cov 0.6
```

### 2.2. Estimate relative abundances
To normalize your data, you will need:

1. Abundance matrix. Should be tab delimited matrix with features (MAGs, genes) in rows and samples in columns. Accepted abundance metrics are: TAD/GEQ, RPKM or reads.
2. MAG metadata. MAG IDs and length (bp) should be provided in 2 columns named as "MAG" and "Length". Genes can be used instead of MAGs but the column header should be kept as "MAG".
3. Metagenome metadata. Metagenome IDs and genome equivalents (GEQ) should be provided in 2 columns named as "Metagenome" and "GEQ". This is only required if data is TAD/GEQ, ignore it for RPKM or reads.
4. .npo files generated after running nonpareil (stand-alone) with each of your metagenomes of interest. These files should be placed on the same folder where you are running the analysis. For instructions about how to run nonpareil visit https://nonpareil.readthedocs.io/en/latest/index.html



Once you have prepared your tables, you are ready to normalize your data in a matter of seconds!


For TAD/GEQ data:
```
Rscript Npc_normalizer.R -a Matrix_TAD_GEQ.txt -f Feature_metadata.txt -m Metagenome_metadata.txt -i TAD_GEQ -n 0.6 
```

For RPKM data:
```
Rscript Npc_normalizer.R -a Matrix_RPKM.txt -f Feature_metadata.txt -i RPKM -n 0.6
```

For reads counts:
```
Rscript Npc_normalizer.R -a Matrix_Reads.txt -f Feature_metadata.txt -i reads -n 0.6
```
