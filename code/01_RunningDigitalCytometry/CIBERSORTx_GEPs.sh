#!/bin/bash

#SBATCH --array=1-2 # there are 2 mixture files

#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_normal
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00

#SBATCH --job-name=cx_GEP
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=elia.alonso@helmholtz-munich.de

# create slurm_out
mkdir -p slurm_out

# Run CIBERSORTx-GEPs module on HPC cluster using SLURM

# 23-09-2025
# Elia Alonso

IMAGE_PATH='/home/idc/elia.alonso/tools/CIBERSORTx/'
USERNAME="elia.alonso@helmholtz-munich.de"
TOKEN=$(cat /home/idc/elia.alonso/tools/CIBERSORTx/mytoken.txt)

INPUT_PATH="/home/idc/elia.alonso/projects/invitro_mixtures/data/input_dig_cyt/CIBERSORTx"
SIG_MATRIX="/home/idc/elia.alonso/projects/invitro_mixtures/data/input_dig_cyt/CIBERSORTx/CIBERSORTx_sig_mat_all_cell_lines.txt"
INPUT_FILES_PATTERN="cx_mixtures_I_"

OUTPUT_PATH="/home/idc/elia.alonso/projects/invitro_mixtures/data/output_dig_cyt/CIBERSORTx"

INPUT_FILES=($(find ${INPUT_PATH} -type f | grep ${INPUT_FILES_PATTERN}))
TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Extract the mixture input file for the current task
MIXTURE=${INPUT_FILES[$TASK_INDEX]}
filename=$(basename "$MIXTURE" .txt)
prefix="GEP"

# Create output folder for this mixture file
OUTPUT_PATH="$OUTPUT_PATH/GEP_output_$filename"
mkdir -p ${OUTPUT_PATH}

# Print task information
echo "Running task $SLURM_ARRAY_TASK_ID with input file: $filename"
echo "Output path: $OUTPUT_PATH"

echo $(date '+%d-%m-%Y')

# go to docker images directory
cd ${IMAGE_PATH}

# start timer
SECONDS=0

# Run CIBERSORTx
apptainer exec -B ${INPUT_PATH}:/src/data -B ${OUTPUT_PATH}:/src/outdir gep_latest.sif Rscript /src/R_modules/CIBERSORTxGEP.R \
    --username ${USERNAME} \
    --token ${TOKEN} \
    --sigmatrix ${SIG_MATRIX} \
    --mixture ${MIXTURE} 

# stop timer
duration=$SECONDS

echo " - Execution time = $((duration / 60)) min $((duration % 60)) sec."
echo " - GEP results: ${OUTPUT_PATH}"
