#! /bin/bash
#
# Script to merge output ROOT files

JOB_ID=26827119

FILE_DIR="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/$JOB_ID"
FILES00=$( find "$FILE_DIR" -path '*000*/00*/*' -name "*.root" )
FILES01=$( find "$FILE_DIR" -path '*000*/01*/*' -name "*.root" )
FILES02=$( find "$FILE_DIR" -path '*000*/02*/*' -name "*.root" )
FILES03=$( find "$FILE_DIR" -path '*000*/03*/*' -name "*.root" )
FILES04=$( find "$FILE_DIR" -path '*000*/04*/*' -name "*.root" )
FILES05=$( find "$FILE_DIR" -path '*000*/05*/*' -name "*.root" )
FILES06=$( find "$FILE_DIR" -path '*000*/06*/*' -name "*.root" )
FILES07=$( find "$FILE_DIR" -path '*000*/07*/*' -name "*.root" )
FILES08=$( find "$FILE_DIR" -path '*000*/08*/*' -name "*.root" )

OUTPUT_DIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/$JOB_ID
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal0.root $FILES00
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal1.root $FILES01
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal2.root $FILES02
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal3.root $FILES03
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal4.root $FILES04
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal5.root $FILES05
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal6.root $FILES06
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal7.root $FILES07
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal8.root $FILES08
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $OUTPUT_DIR/AnalysisResultsFinal*.root
