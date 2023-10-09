#!/bin/bash

# The genes annotated by ANNOVAR sometimes contain \x3b
# Replace with semicolon
sed 's/\\x3b/;/g' tables/9_filtered.csv > tables/10_replace_escape.csv
