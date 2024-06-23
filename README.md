# hiv-cov repository

This repo contains the code for the HIV-COV project. It is somewhat of a mess and still under construction, but should be usable.

## Workflow

1) The `data_setup_BEAST.R` file contains the code to read in the metadata and the .fasta file, process them, and output .fasta files split by lineage with .txt files of dates formatted for BEAUTI2.

2) The `*_beast2` folders contain the BEAST2 .xml files for each lineage. After running the above, they should also contain the split fasta files and files containing dates for BEAUTI2. All BEAST2 results should be saved in these folders to be read in later.

3) The `Analyses_by_lineage` folder contains files to read in the BEAST2 results from the folders in (2), subsample 100 BEAST2 trees, and run TransPhylo on each. This is best to do on the computing cluster. The results are set to be saved to .rds files in a `results` folder. In the same file, they are read back in and probability of individuals in the sample being infection sources are computed and saved to a folder called `prob_source`.

4) The `exploratory.R` file contains some miscellaneous exploratory code for the metadata and the TransPhylo results.
