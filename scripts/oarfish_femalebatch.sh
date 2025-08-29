#!/bin/bash -l
# Batch script to run isoquant.py for hippocampus samples as an array job using a parameter file.
# Request 10 minutes of wallclock time (hh:mm:ss)
#$ -l h_rt=60:00:00

# Request 1 gigabyte of RAM.
#$ -l mem=16G

# Request 20 gigabyte of TMPDIR space.
#$ -l tmpfs=20G

# Set up the job array; here tasks 1-9 correspond to our 4 samples. ##
#$ -t 1-9

# Set the job name.
#$ -N oarfish_000119

# Set the working directory (replace <your_UCL_id> with your UCL user ID).
#$ -wd /home/ucbtpli/Scratch/oarfish_C000119

#$ -j y

# Load environment
conda activate Taf1

# Task-specific parameter extraction
number=$SGE_TASK_ID
paramfile="/home/ucbtpli/Scratch/oarfish_C000119/oarfish_all_params.txt"

# Read the fastq and sample label from the param file
fastq_file="$(sed -n ${number}p $paramfile | awk '{print $1}')"
sample_label="$(sed -n ${number}p $paramfile | awk '{print $2}')"

# Set reference and output paths
REFERENCE="/home/ucbtpli/Scratch/canonical_novel_all_Peihang_nTaf1.fa" 
INDEX_OUT="/home/ucbtpli/Scratch/transcripts_3.mmi"
OUTPUT_DIR="/home/ucbtpli/Scratch/oarfish_C000119/${sample_label}.output"

# Run oarfish
oarfish -j 3 \
  --reads $fastq_file \
  --num-bootstraps 30 \
  --reference $REFERENCE \
  --index-out $INDEX_OUT \
  --seq-tech ont-cdna \
  -o $OUTPUT_DIR \
  --filter-group no-filters \
  --score-threshold 0.9999 \
  --model-coverage \
  --write-assignment-probs

conda deactivate