INFILE="/rstorage/alice/data/LHC16qt/776/LHC16q/000265344/0179/AnalysisResults.root"
INPUT_FILE="/rstorage/alice/data/LHC17pq/448/20-06-2020/448_20200619-0610/unmerged/child_1/2029/AnalysisResults.root"
echo $INPUT_FILE
echo $INPUT_FILE | cut -d/ -f5-11
echo $INFILE
echo $INFILE | cut -d/ -f5-9

