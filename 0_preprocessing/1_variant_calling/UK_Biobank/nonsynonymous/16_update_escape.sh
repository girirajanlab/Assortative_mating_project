#!/bin/bash

# Semicolons in the VCF output were written as \x3b
# Change back to semicolons
sed 's/\\x3b/;/g' tables/15_intracohort_filter.csv > tables/16_replace_escape.csv
