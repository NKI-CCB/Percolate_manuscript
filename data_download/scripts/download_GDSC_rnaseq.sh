data_folder="../../data/"
raw_folder="${data_folder}raw/"
response_folder="${data_folder}response/"

mkdir $raw_folder
mkdir $response_folder

wget -P $raw_folder https://cog.sanger.ac.uk/cmp/download/rnaseq_20191101.zip
unzip "${raw_folder}rnaseq_20191101.zip" -d $raw_folder
rm "${raw_folder}rnaseq_20191101.zip"

wget -P $data_folder https://cog.sanger.ac.uk/cmp/download/model_list_20191104.csv

python ./dependencies/process_GDSC_rnaseq.py

wget -P $response_folder ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_24Jul22.xlsx
wget -P $response_folder ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_24Jul22.xlsx

rm -r $raw_folder