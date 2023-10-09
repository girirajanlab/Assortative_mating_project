#!/bin/bash

# Check if a "large fraction" of samples have a large waviness factor (WF)
# See documentation: https://penncnv.openbioinformatics.org/en/latest/user-guide/test/#gcmodel-adjustment

# Save WF values to file
grep 'WF=' logs/3_autosome_cnvs.log | cut -f 13 -d ' ' > tmp

# ChrX calls do not have a WF

# Filter the results in python
python - << EOF

import pandas as pd

df = pd.read_csv('tmp', sep = '=', names = ['WF', 'value'])
print(df.shape)
df = df[(df.value > 0.04) | (df.value < -0.04)]
print(df.shape)

EOF

rm tmp

# 7 samples had a |WF| > 0.04
# No GC model adjustment will be applied
