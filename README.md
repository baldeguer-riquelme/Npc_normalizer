# Nonpareil coverage normalizer
Repository with scripts to normalize relative abundance (TAD/GEQ, RPKM or recruited reads) to a given nonpareil coverage.

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
cp scripts/Npc_normalizer.R examples
cd examples
Rscript Npc_normalizer.R -a Matrix_TAD_GEQ.txt -f Feature_metadata.txt -m Metagenome_metadata.txt -i TAD_GEQ -n 0.6 
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

### How to run
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
