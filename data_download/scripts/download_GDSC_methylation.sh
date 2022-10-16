data_folder="../../data/"
raw_folder="${data_folder}raw/"

mkdir $raw_folder
wget -P $raw_folder https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/preprocessed/methylation/METH_CELL_DATA.txt.zip
wget -P $raw_folder https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/BEMs/CellLines/CellLines_METH_BEMs.zip
wget -P $raw_folder https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/otherAnnotations/methSampleId_2_cosmicIds.xlsx
unzip "${raw_folder}METH_CELL_DATA.txt.zip" -d $raw_folder

python ./dependencies/process_GDSC_methylation.py
rm -r $raw_folder