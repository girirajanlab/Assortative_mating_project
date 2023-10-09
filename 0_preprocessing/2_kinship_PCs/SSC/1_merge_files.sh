#!/bin/bash

#Merge BED files of all batches together, to create one fileset with all samples
mkdir -p plink_files/1_merge

/data5/software/plink --merge-list list_files/plink_filenames.list --out plink_files/1_merge/SSC_merge
# All SNPs merged - no strand flipping needed
