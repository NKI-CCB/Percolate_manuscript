data_folder="../../data/"
raw_folder="${data_folder}raw/"


mkdir $raw_folder
wget -P $raw_folder https://cog.sanger.ac.uk/cmp/download/mutations_20191101.zip
unzip "${raw_folder}mutations_20191101.zip" -d $raw_folder

python ./dependencies/process_GDSC_mutations.py
rm -r $raw_folder