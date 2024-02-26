#!/bin/bash

# run H-MAGMA on all cell types and the three different traits

celltypes=("Astro_FB" "Astro_PP" "Endothelial" "Exc_L2-3" "Exc_L3-5" "Exc_L4-6_1" "Exc_L4-6_2" "Exc_L4-6_3" "Exc_L5-6" "Exc_L5-6_HTR2C" "In_LAMP5" "In_PVALB_Ba" "In_PVALB_Ch" "In_RELN" "In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")
traits=("SCZ2022" "BIP2021" "MDD")
#ct=("Oligodendrocyte")
#traits=("SCZ2022")
basedir="/Users/nathalie_gerstner/Documents/PostmortemBrain/workspace/scripts/RNA/08_DEanalysis/DESeq2_RELN/"

for trait in "${traits[@]}"; do
  for ct in "${celltypes[@]}"; do
    ${basedir}tables/external/H-MAGMA/magma_v1.10_mac/magma --gene-results ${basedir}tables/external/H-MAGMA/Adult_brain_${trait}.genes.raw --gene-covar ${basedir}tables/GWASenrichment/MAGMAinput/status/${ct}_gene_covar.txt --out ${basedir}tables/GWASenrichment/MAGMAoutput/status/${trait}_${ct}_magmaEnrichment
  done
done

# gene covar
# external/H-MAGMA/magma_v1.10_mac/magma --gene-results external/H-MAGMA/Adult_brain_SCZ2022.genes.raw --gene-covar GWASenrichment/MAGMAinput/Oligodendrocyte_gene_covar.txt --out test_Oligo
