#!/bin/bash

#SBATCH --cpus-per-task=5
#SBATCH --mem=50G
#SBATCH --array=0-21%5
#SBATCH --time=11-00:00:00
#SBATCH --output=slurm_09_7_%A_%a.out
#SBATCH --error=slurm_09_7_%A_%a.err

source ~/.bashrc

conda activate sc-atac

# Run script on all celltypes
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22")

#for i in "${chromosomes[@]}"; do
Rscript /psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/scripts/ATAC/09_7_GeneScoreMatrixToCSV.R "${chromosomes[$SLURM_ARRAY_TASK_ID]}"
#done

conda deactivate