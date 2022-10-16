# Set up the environment
# conda create -y -n data_download python=3.9
eval "$(conda shell.bash hook)"
conda activate data_download
pip install pandas numpy seaborn intervaltree openpyxl

# Create directory
mkdir -p ../../data
for data_folder in count fpkm methylation mutations ns_mutations copy_number_binary copy_number
do
	for filtering in all protein_coding mini_cancer
	do
		mkdir -p ../../data/${filtering}
		mkdir -p ../../data/${filtering}/${data_folder}
	done
done

# Launch download
sh download_GDSC_rnaseq.sh
sh download_GDSC_mutations.sh
sh download_GDSC_copy_number.sh
sh download_GDSC_methylation.sh