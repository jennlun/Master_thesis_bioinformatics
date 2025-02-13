#!/bin/bash
#SBATCH --account=nn8014k
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=normal
#SBATCH --job-name=staralign
#SBATCH --array=0-47

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

module load STAR/2.7.11a-GCC-12.3.0

echo loaded modules

module list

NAMES=($(cat /cluster/work/users/jennlun/rnaseq_sample_names/rnaseq_sample_names_short.list))

echo running STAR alignment on ${NAMES[${SLURM_ARRAY_TASK_ID}]}

cd /cluster/work/users/jennlun/3_map_genome_star/

mkdir ${NAMES[${SLURM_ARRAY_TASK_ID}]}
cd ${NAMES[${SLURM_ARRAY_TASK_ID}]}

/cluster/software/STAR/2.7.11a-GCC-12.3.0/bin/STAR --runThreadN 8 \
 --twopassMode Basic --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate \
 --limitBAMsortRAM 3000000000 --outBAMsortingThreadN 5 --outSAMattributes All \
 --genomeDir /cluster/projects/nn8014k/jenny/3_map_genome_star/ccar_149bp_gtf_index \
 --readFilesIn /cluster/work/users/jennlun/rnaseq_trimmed/${NAMES[${SLURM_ARRAY_TASK_ID}]}_R1_001_val_1.fq.gz /cluster/work/users/jennlun/rnaseq_trimmed/${NAMES[${SLURM_ARRAY_TASK_ID}]}_R2_001_val_2.fq.gz  \
 --readFilesCommand gunzip -c \
 --outFileNamePrefix ${NAMES[${SLURM_ARRAY_TASK_ID}]}_star_map2ccar > ${NAMES[${SLURM_ARRAY_TASK_ID}]}_star_map2ccar.log 2>&1



echo finished STAR alignment