#!/bin/bash -l
# Batch script to run isoquant.py for hippocampus samples as an array job using a parameter file.
# Request 10 minutes of wallclock time (hh:mm:ss)
#$ -l h_rt=60:00:0

# Request 1 gigabyte of RAM.
#$ -l mem=16G

# Request 20 gigabyte of TMPDIR space.
#$ -l tmpfs=20G

# Set up the job array; here tasks 1-9 correspond to our 4 samples.
#$ -t 1-9

# Set the job name.
#$ -N C57Taf1_6female_isoquant

# Set the working directory (replace <your_UCL_id> with your UCL user ID).
#$ -wd /home/ucbtpli/Scratch/C000119_femalemouse_isoquant

#$ -j y

conda activate Taf1

# Get the current task number from the array.
number=$SGE_TASK_ID
paramfile=/home/ucbtpli/Scratch/C000119_femalemouse_isoquant/paras_female3.txt 
# Extract parameters from the parameter file.
sample_path="$(sed -n ${number}p $paramfile | awk '{print $1}')"
sample_label="$(sed -n ${number}p $paramfile | awk '{print $2}')"

# Define common variables.
REF="/home/ucbtpli/Scratch/C000119_femalemouse_isoquant/GRCm39.primary_assembly.genome.fa"
GTF="/home/ucbtpli/Scratch/C000119_femalemouse_isoquant/gencode.vM36.annotation.db"
base_output_dir="/home/ucbtpli/Scratch/C000119_femalemouse_isoquant" 
output_dir="${base_output_dir}/${sample_label}"

echo "Processing sample: ${sample_label}"

isoquant.py --reference $REF \
--bam $sample_path \
--genedb $GTF \
--sqanti_output \
--data_type nanopore -o $output_dir

conda deactivate
