#! /bin/bash

INITDIR=$(pwd)
OUTDIR=testmc/embedding/ #CHANGE!!!
mkdir $OUTDIR

# UPDATE CONFIG FILE IF NEEDED!!!

for i in {2..2} #7
do
    [ ! -d "$OUTDIR/$i" ] && mkdir $OUTDIR/$i
    # MCTESTFILE="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/data/LHC18b8/569/LHC18b8_fast/$i/282008/0001/AnalysisResults.root"
    # MCTESTFILE="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/data/LHC18b8/569/LHC18b8_fast/19/28230$i/0001/AnalysisResults.root"
    MCTESTFILE="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/data/LHC18b8/569/LHC18b8_cent_woSDD/8/282341/0007/AnalysisResults.root"
    MCTESTDIR=$OUTDIR/$i
    echo $MCTESTDIR
    python process_mc_ENC.py -c $AAPYJ/config/ENC/pPb/process_pPb.yaml -f $MCTESTFILE -o $MCTESTDIR
done

# cd $OUTDIR
# rm AnalysisResultsFinal.root
# #python $AAPYJ/slurm/utils/james/scaleHistograms.py -c $AAPYJ/slurm/utils/anjali/scaleFactors.yaml
# hadd -f AnalysisResultsFinal.root */*.root
# cd $INITDIR
