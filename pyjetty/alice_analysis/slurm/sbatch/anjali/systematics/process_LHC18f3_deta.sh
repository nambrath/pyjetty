#! /bin/bash

# This script takes an input file path as an argument, and runs a python script to 
# process the input file and write an output ROOT file.
# The main use is to give this script to a slurm script.

# Take two command line arguments -- (1) input file path, (2) output dir prefix
if [ "$1" != "" ]; then
  INPUT_FILE=$1
  #echo "Input file: $INPUT_FILE"
else
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  JOB_ID=$2
  echo "Job ID: $JOB_ID"
else 
  echo "Wrong command line arguments"
fi

if [ "$3" != "" ]; then
  TASK_ID=$3
  echo "Task ID: $TASK_ID"
else
  echo "Wrong command line arguments"
fi

# Define output path from relevant sub-path of input file
OUTPUT_PREFIX="AnalysisResults/anjali/$JOB_ID"
# Note: suffix depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f7-11)
echo "Output suffix: $OUTPUT_SUFFIX"
#echo $OUTPUT_SUFFIX
OUTPUT_DIR="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/$OUTPUT_PREFIX/$OUTPUT_SUFFIX"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
module use /global/cfs/cdirs/alice/heppy_soft/yasp/software/modules
module load cmake gsl root/6.28.00 HepMC2/2.06.11 LHAPDF6/6.5.3 pcre2/default swig/4.1.1 HepMC3/3.2.5
module use $workdir/pyjetty/modules
module load pyjetty

# Run python script via pipenv
cd /global/cfs/cdirs/alice/anjali/mypyjetty/pyjetty/pyjetty/alice_analysis
python process/user/anjali/process_mc_ENC.py -c config/ENC/pPb/systematics/process_pPb_mc_deta.yaml -f $INPUT_FILE -o $OUTPUT_DIR

# Move stdout to appropriate folder
mv /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/slurm-${JOB_ID}_${TASK_ID}.out /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/${JOB_ID}/
