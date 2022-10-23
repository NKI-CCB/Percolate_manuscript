#!/bin/sh
# LAUNCH ESTIMATION OF COMPONENTS USING AIC ON A GRID SEARCH

script_folder=./model_training/

# Results folder
output_folder=./output/model_selection/
mkdir -p ${output_folder}

for data_type in MUT CNA METH #GE
do
	echo "START ${data_type}"
	mkdir -p ${output_folder}/${data_type} ${output_folder}/${data_type}/tmp/
	python -Xfaulthandler ${script_folder}/divide_data.py -o ${output_folder}/${data_type} -d ${data_type}
	python ${script_folder}/compute_likelihood_gridsearch.py \
		-o ${output_folder}/${data_type} -t ${output_folder}/${data_type}/tmp/ -d ${data_type}
done
