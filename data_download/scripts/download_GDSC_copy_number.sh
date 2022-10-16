data_folder="../../data/"
raw_folder="${data_folder}raw/"

mkdir $raw_folder

wget -P $raw_folder https://cog.sanger.ac.uk/cmp/download/cnv_20191101.zip
unzip "${raw_folder}cnv_20191101.zip" -d $raw_folder

python ./dependencies/process_GDSC_copy_number.py
rm -r $raw_folder