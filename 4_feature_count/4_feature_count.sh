#!/bin/bash
#SBATCH --account=nn8014k
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=featureCounts


set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load Subread/2.0.4-GCC-11.3.0 #load module

echo running featureCounts

featureCounts -T 5 -O -C -p -s 2 --countReadPairs -t exon -g gene_id -a /cluster/projects/nn8014k/jenny/3_map_genome_star/ccar_149bp_gtf_index -o rnaseq_map2ccar_featureCounts.txt *_star_map2ccarAligned.sortedByCoord.out.bam


echo finished featureCounts
