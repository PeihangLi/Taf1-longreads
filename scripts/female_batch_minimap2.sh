#!/bin/bash

# set file paths
GENOME="/Users/lipeihang/Desktop/nanopore_aligned_data/GRCm39.primary_assembly.genome.fa"
ANNOTATION="/Users/lipeihang/Desktop/nanopore_aligned_data/gencode.vM36.annotation.bed12"
MICRO_EXON_BED="/Users/lipeihang/Desktop/nanopore_aligned_data/micro_exon_34.bed"

# define folders to process
FOLDERS=("912979_Hippocampus" "912979_Striatum" "912985_Cerebellum" "912985_Cortex" "912985_Hippoacmpus" "912985_Striatum" "912998_Striatum")

# Loop through all .fastq files
for folder in "${FOLDERS[@]}"; do
    # Ensure folder exists
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"

        # Find all .fastq files in current folder
        for fastq_file in "$folder"/*.fastq; do
            # Ensure file exists (prevent error if directory is empty)
            [ -e "$fastq_file" ] || continue
            
            # Get sample name (remove path and .fastq extension)
            sample_name=$(basename "$fastq_file" .fastq)
            
            echo "Processing sample: $sample_name"

            # alignment using minimap2
            minimap2 -ax splice "$GENOME" "$fastq_file" > "$folder/${sample_name}.sam"

            # convert SAM to BAM
            samtools view -S -b "$folder/${sample_name}.sam" > "$folder/${sample_name}.bam"

            # run MisER for micro-exon alignment
            MisER -c4 -s10 -d0.2 -f20 --allTranscripts --fixFlankLen --debugMode \
                  "$folder/${sample_name}.bam" "$GENOME" "$ANNOTATION" "$MICRO_EXON_BED" \
                  -o "$folder/${sample_name}_M34.bam"

            # sort BAM
            samtools sort "$folder/${sample_name}_M34.bam" -o "$folder/${sample_name}_M34sorted.bam"

            # index BAM
            samtools index "$folder/${sample_name}_M34sorted.bam"

            echo "Finished processing: $sample_name"
        done
    else
        echo "Warning: Folder $folder does not exist!"
    fi
done

echo "All samples processed successfully!"

