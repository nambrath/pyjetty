#! /bin/bash

# Script to merge output ROOT files
JOB_ID=26090830
OUTPUT_DIR="$MYAR/$JOB_ID"

# command line arguments
if [ "$1" != "" ]; then
  MERGE_JOB_ID=$1
  echo "Merge Job ID: $MERGE_JOB_ID"
else
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  BIN=$2
  echo "Bin: $BIN"
else
  echo "Wrong command line arguments"
fi

# Load modules
module use /global/cfs/cdirs/alice/heppy_soft/yasp/software/modules
module load cmake gsl root/6.28.00 HepMC2/2.06.11 LHAPDF6/6.5.3 pcre2/default swig/4.1.1 HepMC3/3.2.5
module use $workdir/pyjetty/modules
module load pyjetty

# Merge all output files from each pt-hat bin
FILE_DIR_BASE=$MYAR/$JOB_ID
FILES=$( find ${FILE_DIR_BASE}/*/${BIN}/*/* -name "*.root" )

OUT_DIR_BASE=$MYAR/$JOB_ID
mkdir -p ${OUT_DIR_BASE}/Stage0/${BIN}
hadd -f -j 10 ${OUT_DIR_BASE}/Stage0/${BIN}/AnalysisResults.root $FILES

# Move stdout to appropriate folder
mv $MYAR/slurm-${MERGE_JOB_ID}_${BIN}.out $MYAR/${JOB_ID}/
