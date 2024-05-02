# Nonpareil coverage normalizer
Repository with scripts to normalize relative abundance (TAD/GEQ, RPKM or recruited reads) at a given nonpareil coverage

### Install dependencies
```
conda create -n npc_norm -c conda-forge -c bioconda r-tidyverse r-reshape2 nonpareil r-roxygen2
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

### How to run
To normalize your data, you will need:

1. Abundance matrix. Should be tab delimited matrix with features (MAGs, genes) in rows and samples in columns. Accepted abundance metrics are: TAD/GEQ, RPKM or reads.
2. MAG metadata. MAG IDs and length (bp) should be provided in 2 columns named as "MAG" and "Length". Genes can be used instead of MAGs but the column header should be kept as "MAG".
3. Metagenome metadata. Metagenome IDs and genome equivalents (GEQ) should be provided in 2 columns named as "Metagenome" and "GEQ". This is only required if data is TAD/GEQ, ignore it for RPKM or reads.
4. .npo files generated after running nonpareil (stand-alone) with each of your metagenomes of interest. These files should be placed on the same folder where you are running the analysis. For instructions about how to run nonpareil visit https://nonpareil.readthedocs.io/en/latest/index.html



Once you have prepared your tables, you are ready to normalize your data in a matter of seconds!


For TAD/GEQ data:
```
Rscript Npc_normalizer.R
```

For RPKM data:
```
Rscript Npc_normalizer.R
```

For reads counts:
```
Rscript Npc_normalizer.R
```
