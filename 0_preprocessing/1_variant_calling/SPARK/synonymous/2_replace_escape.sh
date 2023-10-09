#!/bin/bash

# The genes annotated by ANNOVAR sometimes contain \x3b
# Replace with semicolon
sed 's/\\x3b/;/g' tables/1_synonymous.csv > tables/2_replace_escape.csv
