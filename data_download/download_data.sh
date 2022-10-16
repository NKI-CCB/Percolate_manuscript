# Set up the environment
# conda create -y -n data_download python=3.9
eval "$(conda shell.bash hook)"
conda activate data_download
# pip install pandas numpy seaborn

sh download_GDSC.sh