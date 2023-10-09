#!/bin/bash

#Merge BED files of Data1 and Data2 batches together, to create one fileset with all samples
mkdir -p plink_files/1_merge

/data5/software/plink --merge-list list_files/plink_filenames.list --out plink_files/1_merge/SPARK_data12
# Perform strand-flipping
/data5/software/plink --bfile array/genotype/SPARK.iWES_v1.array.2022_02 --flip plink_files/1_merge/SPARK_data12.missnp --make-bed --out plink_files/1_merge/data1-2_flip_test
/data5/software/plink --bfile array/genotype/SPARK.iWES_v1.array.2022_02 --bmerge plink_files/1_merge/data1-2_flip_test --make-bed --out plink_files/1_merge/merged_test
#Exclude variants that fail strand flipping
/data5/software/plink --bfile array/genotype/SPARK.iWES_v1.array.2022_02 --exclude plink_files/1_merge/merged_test-merge.missnp --make-bed --out plink_files/1_merge/SPARK_data2_temp
/data5/software/plink --bfile plink_files/1_merge/data1-2_flip_test --exclude plink_files/1_merge/merged_test-merge.missnp --make-bed --out plink_files/1_merge/SPARK_data1_temp
#Perform final merge
/data5/software/plink --bfile plink_files/1_merge/SPARK_data2_temp --bmerge plink_files/1_merge/SPARK_data1_temp --make-bed --out plink_files/1_merge/SPARK_merged
