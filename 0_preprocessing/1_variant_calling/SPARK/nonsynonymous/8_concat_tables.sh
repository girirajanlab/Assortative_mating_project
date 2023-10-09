#!/bin/bash

# Combine all the tables into a single file
cat tables/7_filter_variants/* >> tables/8_combined.tsv
