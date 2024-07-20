#! /bin/bash

# INITDIR=$(pwd)
# OUTDIR=~/embedROTtest #CHANGE!!!
# mkdir $OUTDIR
# rm -r $OUTDIR/*

rm ~/FILELIST18b8.txt 

for i in {1..20}
do
    [ ! -d "$OUTDIR/$i" ] && mkdir $OUTDIR/$i
    find /rstorage/alice/data/LHC18b8/569/LHC18b8_fast/$i -name "AnalysisResults.root" > ~/testfile.txt
    a=($(wc ~/testfile.txt))
    echo $i $a
    let filenum=${a[0]}/10
    head -$filenum ~/testfile.txt >> ~/FILELIST18b8.txt
done



# cd $OUTDIR
# rm AnalysisResultsFinal.root
# python $AAPYJ/slurm/utils/james/scaleHistograms.py -c /rstorage/alice/data/LHC18b8/scaleFactors.yaml
# hadd -f AnalysisResultsFinal.root */*/*.root
# cd $INITDIR
