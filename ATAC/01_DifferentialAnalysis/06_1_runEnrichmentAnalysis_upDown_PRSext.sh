#!/bin/bash

#SBATCH --partition=p.psycl
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time=11-00:00:00
#SBATCH --output=slurm_06_1_%j.out
#SBATCH --error=slurm_06_1_%j.err

source ~/.bashrc

conda activate clusterProfiler

prstrait=("crossDisorder2019" "BIP2021" "MDD" "SCZ2022" "height")

for p in "${prstrait[@]}"; do
    Rscript /psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/scripts/ATAC/02_DifferentialAnalysis/06_1_EnrichmentAnalysis_upDown_PRSext.R "$p"
done

conda deactivate