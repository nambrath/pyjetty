# DATATESTFILE="/rstorage/alice/data/LHC16qt/776/LHC16q_cent_woSDD/000265344/0179/AnalysisResults.root"
# DATATEST_DIR="~/testdata"

# python process_data_ENC.py -c $AAPYJ/config/ENC/pPb/process_pPb.yaml -f $DATATESTFILE -o $DATATEST_DIR



# run on a few 
INITDIR=$(pwd)
DATATEST_DIR=testdata/ppbperp
mkdir $DATATEST_DIR
FILE_PATHS='/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/data/LHC16qt/776/files.txt'
# for (( i = 1; i <= 20; i++ ))
for (( i = 1; i <= 1; i++ ))
do
  FILE=/global/cfs/projectdirs/alice/alicepro/hiccup/$(sed -n "$i"p $FILE_PATHS)
  [ ! -d "$DATATEST_DIR/$i" ] && mkdir $DATATEST_DIR/$i
  FILEDIR=$DATATEST_DIR/$i
  echo $FILEDIR
  python process_data_ENC.py -c $AAPYJ/config/ENC/pPb/process_pPb.yaml -f $FILE -o $FILEDIR
  # python process_data_ENC.py -c $AAPYJ/config/ENC/pPb/systematics/process_pPb_deta.yaml -f $FILE -o $FILEDIR
done

cd $DATATEST_DIR
rm AnalysisResultsFinal.root
hadd -f AnalysisResultsFinal.root */*.root
cd $INITDIR
