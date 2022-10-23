# Percolate_manuscript

This GitHub repository contains the scripts used to produce the results for the manuscript <b><a href="https://www.biorxiv.org/content/10.1101/2022.09.11.507473v1">Designing DNA-based predictors of drug response using the signal joint with gene expression</a></b>.
<br/>
The scripts require anaconda (or mini-conda).

## Conda environment

Using <b>create_environment.sh</b> would create the <em>percolate_manuscript</em> environment with all required packages installed.

## Reproducing results
To reproduce the results presented in the manuscript, you can follow these steps.

### First step: Downloading data

Using <b>data_download/scripts/download_GDSC.sh</b> will automatically download and process all the data needed for reproducing the different figures. Downloaded and processed files will appear in the <b>data</b> folder.

## Second step: model selection and training of Percolate models.

Using <b>sh model_training/launch_GDSC_estimation_components_gridsearchAIC.sh</b> would launch the model selection by Grid Search (AIC), train the different GLM-PCA models and align the models by Percolate. Results are saved in <b>output</b>.

## Citation

If you use scripts figuring in this repo, please cite <a href="https://www.biorxiv.org/content/10.1101/2022.09.11.507473v1"><b>Designing DNA-based predictors of drug response using the signal joint with gene expression</b>, Mourragui et al 2022, Biorxiv.</a>