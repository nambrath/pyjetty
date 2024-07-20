#! /bin/bash
#
# Script to merge output ROOT files from all pt-hat bins together, in stages
JOB_ID=26090830
FILE_DIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/$JOB_ID
OUTPUT_DIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/$JOB_ID

# Merge all output files from each pt-hat bin
mkdir -p $OUTPUT_DIR
hadd -f -j 10 $OUTPUT_DIR/AnalysisResultsFinal_scaled.root $FILE_DIR/Stage0/*/*.root
