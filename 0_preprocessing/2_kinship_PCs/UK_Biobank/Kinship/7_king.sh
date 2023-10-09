#!/bin/bash

#Use KING software to calculate pairwise kinship coefficients for all parents (http://people.virginia.edu/~wc9c/KING/manual.html)
/data5/software/king -b plink_files/6_final_filter/UKB_king_input.bed --kinship --degree 1

