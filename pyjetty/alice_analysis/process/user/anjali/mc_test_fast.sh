#! /bin/bash

INITDIR=$(pwd)
OUTDIR=testmc/dpmjet/JETPTMIG #CHANGE!!!
mkdir $OUTDIR

# UPDATE CONFIG FILE IF NEEDED!!!

# MCTESTFILE="/global/cfs/cdirs/alice/anjali/LHC18f3/794/cent_woSDD_1/265500/0037/AnalysisResults.root"
# python process_mc_ENC.py -c $AAPYJ/config/ENC/pPb/process_pPb_mc.yaml -f $MCTESTFILE -o $OUTDIR


for i in {01..01} #20}
do
    [ ! -d "$OUTDIR/$i" ] && mkdir $OUTDIR/$i
    # MCTESTFILE="/global/cfs/cdirs/alice/anjali/LHC20f11c2/793/fast/265383/00$i/AnalysisResults.root"
    MCTESTFILE="/global/cfs/cdirs/alice/anjali/LHC18f3/794/fast_1/265334/00$i/AnalysisResults.root"
    MCTESTDIR=$OUTDIR/$i
    echo $MCTESTDIR
    python process_mc_ENC.py -c $AAPYJ/config/ENC/pPb/process_pPb_mc.yaml -f $MCTESTFILE -o $MCTESTDIR
    # python process_mc_ENC.py -c $AAPYJ/config/ENC/pPb/systematics/process_pPb_mc_norho.yaml -f $MCTESTFILE -o $MCTESTDIR
done

# cd $OUTDIR
# rm AnalysisResultsFinal.root
# hadd -f AnalysisResultsFinal.root */*.root
# cd $INITDIR
