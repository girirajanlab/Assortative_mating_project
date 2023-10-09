#!/bin/bash

# Calculate ROH for samples
# Use same parameters as Gamsiz et al. Am J Hum Genet 2013 (https://pubmed.ncbi.nlm.nih.gov/23830515/)
# Allow 1 heterozygous and 5 missing calls per window
/data5/software/plink --bfile Analysis_files/2_subset_probands/2_roh_samples --homozyg --homozyg-kb 1000 --homozyg-window-het 1 --homozyg-window-missing 5 --out Result_files/3.1_1kb_window
