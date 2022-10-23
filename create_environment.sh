# Set up the environment
conda create -y -n percolate_manuscript python=3.9
eval "$(conda shell.bash hook)"
conda activate percolate_manuscript
pip install pandas numpy seaborn mctorch torch upsetplot matplotlib==3.3.2 rpy2 hyperopt
mamba install -c r r=3.5.1
mamba install -y -c anaconda ipykernel

git clone https://github.com/saroudant/Percolate
cd Percolate
pip install .


python -m ipykernel install --user --name=percolate_manuscript

python install_Rdep.py