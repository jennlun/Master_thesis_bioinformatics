#!/bin/bash

#SBATCH --account=nn8014k
#SBATCH --job-name=TrimGalore_cc
#SBATCH --time=22:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --array=0-95

# load the needed modules
module purge 
ml load Trim_Galore/0.6.10-GCCcore-11.3.0
ml load cutadapt/4.2-GCCcore-11.3.0
ml load FastQC/0.11.9-Java-11

# show me the loaded modules in a list
ml list

# change directory to the output directory in work 
cd /cluster/work/users/jennlun/2_trim_galore/

# Command of the program 
input1_forward=("10-G-N6_S173_L007_R1_001.fastq.gz" \
         "11-G-N7_S174_L007_R1_001.fastq.gz" \
         "12-G-N8_S175_L007_R1_001.fastq.gz" \
         "13-G-R11_S176_L007_R1_001.fastq.gz" \
         "14-G-R12_S177_L007_R1_001.fastq.gz" \
         "15-G-R13_S178_L007_R1_001.fastq.gz" \
         "16-G-R14_S179_L007_R1_001.fastq.gz" \
         "17-G-R17_S180_L007_R1_001.fastq.gz" \
         "18-G-R18_S181_L007_R1_001.fastq.gz" \
         "19-G-R22_S182_L007_R1_001.fastq.gz" \
         "1-G-A1_S165_L007_R1_001.fastq.gz" \
         "20-G-R23_S183_L007_R1_001.fastq.gz" \
         "21-G-R24_S184_L007_R1_001.fastq.gz" \
         "22-G-R25_S185_L007_R1_001.fastq.gz" \
         "23-G-R27_S186_L007_R1_001.fastq.gz" \
         "24-G-R28_S187_L007_R1_001.fastq.gz" \
         "25-L-A1_S188_L007_R1_001.fastq.gz" \
         "26-L-A2_S189_L007_R1_001.fastq.gz" \
         "27-L-A3_S190_L007_R1_001.fastq.gz" \
         "28-L-A4_S191_L007_R1_001.fastq.gz" \
         "29-L-A5_S192_L007_R1_001.fastq.gz" \
         "2-G-A2_S166_L007_R1_001.fastq.gz" \
         "30-L-A6_S193_L007_R1_001.fastq.gz" \
         "31-L-N1_S194_L007_R1_001.fastq.gz" \
         "32-L-N2_S195_L007_R1_001.fastq.gz" \
         "33-L-N3_S196_L007_R1_001.fastq.gz" \
         "34-L-N4_S197_L007_R1_001.fastq.gz" \
         "35-L-N5_S198_L007_R1_001.fastq.gz" \
         "36-L-N6_S199_L007_R1_001.fastq.gz" \
         "37-L-R11_S200_L007_R1_001.fastq.gz" \
         "38-L-R12_S201_L007_R1_001.fastq.gz" \
         "39-L-R13_S202_L007_R1_001.fastq.gz" \
         "3-G-A3_S167_L007_R1_001.fastq.gz" \
         "40-L-R14_S203_L007_R1_001.fastq.gz" \
         "41-L-R15_S204_L007_R1_001.fastq.gz" \
         "42-L-R16_S205_L007_R1_001.fastq.gz" \
         "43-L-R21_S206_L007_R1_001.fastq.gz" \
         "44-L-R22_S207_L007_R1_001.fastq.gz" \
         "45-L-R23_S208_L007_R1_001.fastq.gz" \
         "46-L-R24_S209_L007_R1_001.fastq.gz" \
         "47-L-R25_S210_L007_R1_001.fastq.gz" \
         "48-L-R26_S211_L007_R1_001.fastq.gz" \
         "4-G-A5_S168_L007_R1_001.fastq.gz" \
         "5-G-A7_S169_L007_R1_001.fastq.gz" \
         "6-G-A8_S170_L007_R1_001.fastq.gz" \
         "7-G-N2_S171_L007_R1_001.fastq.gz" \
         "8-G-N3_S212_L007_R1_001.fastq.gz" \
         "9-G-N4_S172_L007_R1_001.fastq.gz")

input1_reverse=("10-G-N6_S173_L007_R2_001.fastq.gz" \
         "11-G-N7_S174_L007_R2_001.fastq.gz" \
         "12-G-N8_S175_L007_R2_001.fastq.gz" \
         "13-G-R11_S176_L007_R2_001.fastq.gz" \
         "14-G-R12_S177_L007_R2_001.fastq.gz" \
         "15-G-R13_S178_L007_R2_001.fastq.gz" \
         "16-G-R14_S179_L007_R2_001.fastq.gz" \
         "17-G-R17_S180_L007_R2_001.fastq.gz" \
         "18-G-R18_S181_L007_R2_001.fastq.gz" \
         "19-G-R22_S182_L007_R2_001.fastq.gz" \
         "1-G-A1_S165_L007_R2_001.fastq.gz" \
         "20-G-R23_S183_L007_R2_001.fastq.gz" \
         "21-G-R24_S184_L007_R2_001.fastq.gz" \
         "22-G-R25_S185_L007_R2_001.fastq.gz" \
         "23-G-R27_S186_L007_R2_001.fastq.gz" \
         "24-G-R28_S187_L007_R2_001.fastq.gz" \
         "25-L-A1_S188_L007_R2_001.fastq.gz" \
         "26-L-A2_S189_L007_R2_001.fastq.gz" \
         "27-L-A3_S190_L007_R2_001.fastq.gz" \
         "28-L-A4_S191_L007_R2_001.fastq.gz" \
         "29-L-A5_S192_L007_R2_001.fastq.gz" \
         "2-G-A2_S166_L007_R2_001.fastq.gz" \
         "30-L-A6_S193_L007_R2_001.fastq.gz" \
         "31-L-N1_S194_L007_R2_001.fastq.gz" \
         "32-L-N2_S195_L007_R2_001.fastq.gz" \
         "33-L-N3_S196_L007_R2_001.fastq.gz" \
         "34-L-N4_S197_L007_R2_001.fastq.gz" \
         "35-L-N5_S198_L007_R2_001.fastq.gz" \
         "36-L-N6_S199_L007_R2_001.fastq.gz" \
         "37-L-R11_S200_L007_R2_001.fastq.gz" \
         "38-L-R12_S201_L007_R2_001.fastq.gz" \
         "39-L-R13_S202_L007_R2_001.fastq.gz" \
         "3-G-A3_S167_L007_R2_001.fastq.gz" \
         "40-L-R14_S203_L007_R2_001.fastq.gz" \
         "41-L-R15_S204_L007_R2_001.fastq.gz" \
         "42-L-R16_S205_L007_R2_001.fastq.gz" \
         "43-L-R21_S206_L007_R2_001.fastq.gz" \
         "44-L-R22_S207_L007_R2_001.fastq.gz" \
         "45-L-R23_S208_L007_R2_001.fastq.gz" \
         "46-L-R24_S209_L007_R2_001.fastq.gz" \
         "47-L-R25_S210_L007_R2_001.fastq.gz" \
         "48-L-R26_S211_L007_R2_001.fastq.gz" \
         "4-G-A5_S168_L007_R2_001.fastq.gz" \
         "5-G-A7_S169_L007_R2_001.fastq.gz" \
         "6-G-A8_S170_L007_R2_001.fastq.gz" \
         "7-G-N2_S171_L007_R2_001.fastq.gz" \
         "8-G-N3_S212_L007_R2_001.fastq.gz" \
         "9-G-N4_S172_L007_R2_001.fastq.gz")

output=("GN6" \
         "GN7" \
         "GN8" \
         "GR11" \
         "GR12" \
         "GR13" \
         "GR14" \
         "GR17" \
         "GR18" \
         "GR22" \
         "GA1" \
         "GR23" \
         "GR24" \
         "GR25" \
         "GR27" \
         "GR28" \
         "LA1" \
         "LA2" \
         "LA3" \
         "LA4" \
         "LA5" \
         "GA2" \
         "LA6" \
         "LN1" \
         "LN2" \
         "LN3" \
         "LN4" \
         "LN5" \
         "LN6" \
         "LR11" \
         "LR12" \
         "LR13" \
         "GA3" \
         "LR14" \
         "LR15" \
         "LR16" \
         "LR21" \
         "LR22" \
         "LR23" \
         "LR24" \
         "LR25" \
         "LR26" \
         "GA5" \
         "GA7" \
         "GA8" \
         "GN2" \
         "GN3" \
         "GN4")    

trim_galore \
--cores 2 \
--paired \
--quality 20 \
--fastqc \
-o /cluster/work/users/jennlun/2_trim_galore/${output[$SLURM_ARRAY_TASK_ID]} \
 /cluster/work/users/jennlun/2024_raw_rnaseq_data/240506_LH00534.B.Project_Lundeberg-RNA1-2024-03-18/${input1_forward[$SLURM_ARRAY_TASK_ID]} \
 /cluster/work/users/jennlun/2024_raw_rnaseq_data/240506_LH00534.B.Project_Lundeberg-RNA1-2024-03-18/${input1_reverse[$SLURM_ARRAY_TASK_ID]}


# to close everything 
ml purge 