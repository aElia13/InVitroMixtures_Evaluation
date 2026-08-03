#!/bin/bash

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16
#SBATCH --time=0:10:00
#SBATCH --ntasks=1

#SBATCH --job-name=cx_Fractions
#SBATCH --output=slurm_out/%x_%j_output.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=elia.alonso@helmholtz-munich.de

# create slum_out
mkdir -p slurm_out

# Run CIBERSORTx-Fractions module to build signature matrix

# 23-09-2025
# Elia Alonso

IMAGE_PATH='/home/idc/elia.alonso/tools/CIBERSORTx/'
USERNAME="elia.alonso@helmholtz-munich.de"
TOKEN=$(cat /home/idc/elia.alonso/tools/CIBERSORTx/mytoken.txt)

INPUT_PATH='/home/idc/elia.alonso/projects/invitro_mixtures/data/input_dig_cyt/CIBERSORTx'
REF_SAMPLE_FILE='/home/idc/elia.alonso/projects/invitro_mixtures/data/input_dig_cyt/CIBERSORTx/cx_reference_sample_file_all_cell_lines.txt'
PHENO_CLASSES_FILE='/home/idc/elia.alonso/projects/invitro_mixtures/data/input_dig_cyt/CIBERSORTx/cx_phenotype_classes_file_all_cell_lines.txt'

OUTPUT_PATH="/home/idc/elia.alonso/projects/invitro_mixtures/data/output_dig_cyt/CIBERSORTx/Fractions_output_all_cell_lines"
mkdir -p ${OUTPUT_PATH}

echo $(date '+%d-%m-%Y')

# go to docker images directory
cd ${IMAGE_PATH}

# Run CIBERSORTx
apptainer exec -B ${INPUT_PATH}:/src/data -B ${OUTPUT_PATH}:/src/outdir cibersortx_fractions.sif /src/CIBERSORTxFractions \
--username ${USERNAME} \
--token ${TOKEN} \
--single_cell FALSE \
--refsample ${REF_SAMPLE_FILE} \
--phenoclasses ${PHENO_CLASSES_FILE} \
--G.min 50 \
--G.max 150 \
--q.value 0.3 \
--verbose TRUE 

echo " - Fractions results: ${OUTPUT_PATH}"
