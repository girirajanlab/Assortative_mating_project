#!/bin/bash
OUT_DIR=vcfs/2_quality_filtered/

for i in `seq 0 9`
do
  mkdir -p $OUT_DIR/$i
done
