#!/bin/bash

#SBATCH --partition=p.psycl
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --time=11-00:00:00
#SBATCH --output=slurm_03_2_%j.out
#SBATCH --error=slurm_03_2_%j.err

source ~/.bashrc

conda activate sc-atac

# Run script on all celltypes
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22")

for i in "${chromosomes[@]}"; do
  Rscript /psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/scripts/ATAC/02_DifferentialAnalysis/03_2_GeneratePseudobulkGeneScoreMatrix.R "$i"
done

conda deactivate