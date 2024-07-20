#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=26827125

FILE_DIR="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/$JOB_ID"
FILES=$( find "$FILE_DIR" -name "*.root" )

OUTPUT_DIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/$JOB_ID
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
