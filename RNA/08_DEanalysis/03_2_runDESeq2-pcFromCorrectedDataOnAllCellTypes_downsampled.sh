ct=("Astro_FB" "Astro_PP" "Endothelial" "Exc_L2-3" "Exc_L3-5" "Exc_L4-6_1" "Exc_L4-6_2" "Exc_L4-6_3" "Exc_L5-6" "Exc_L5-6_HTR2C" "In_LAMP5" "In_PVALB_Ba" "In_PVALB_Ch" "In_RELN" "In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")
#ct=("In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")
perc=(25 50)
# perc=(75)

for p in "${perc[@]}"; do
  for i in "${ct[@]}"; do
  Rscript ~/Documents/PostmortemBrain/workspace/scripts/RNA/08_DEanalysis/DESeq2_RELN/03_2_DESeq2-pcFromCorrectedDataOnAllCellTypes_downsampled.R "$i" "$p"
  done
done