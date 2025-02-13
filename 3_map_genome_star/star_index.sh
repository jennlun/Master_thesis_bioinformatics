#!/bin/bash
#SBATCH --account=nn8014k
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --partition=bigmem
#SBATCH --job-name=starindex

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

module load STAR/2.7.11a-GCC-12.3.0

echo Loaded modules

module list

echo generating STAR genome index

/cluster/software/STAR/2.7.11a-GCC-12.3.0/bin/STAR --runMode genomeGenerate --genomeDir ccar_149bp_gtf_index --genomeFastaFiles /cluster/projects/nn8014k/ccar_latest_annotation/ccar_genome_v1_262scaffolds.fasta --sjdbGTFfile /cluster/projects/nn8014k/ccar_latest_annotation/carcar_annotation_v5.gtf --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 149 > ccar_star_index149.log 2>&1

echo finished generating STAR genome index
