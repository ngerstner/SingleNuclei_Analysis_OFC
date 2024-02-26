#!/bin/bash

#SBATCH --partition=p.psycl
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --time=11-00:00:00
#SBATCH --output=slurm_03_3e_%j.out
#SBATCH --error=slurm_03_3e_%j.err

source ~/.bashrc

conda activate DESeq2

# Run script on all celltypes
ct=("Astro_FB" "Astro_PP" "Endothelial" "Exc_L2-3" "Exc_L3-5" "Exc_L4-6_1" "Exc_L4-6_2" "In_PVALB_Ba" "In_PVALB_Ch" "In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")

for i in "${ct[@]}"; do
  #Rscript /psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/scripts/ATAC/02_DifferentialAnalysis/03_1_DESeq2-pcFromCorrectedDataOnAllCellTypes_newFiltering.R "$i"
  Rscript /psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/scripts/ATAC/02_DifferentialAnalysis/03_3e_CorrelationRNAandPeak.R "$i"
done

conda deactivate
