ct=("Astro_FB" "Astro_PP" "Endothelial" "Exc_L2-3" "Exc_L3-5" "Exc_L4-6_1" "Exc_L4-6_2" "Exc_L4-6_3" "Exc_L5-6" "Exc_L5-6_HTR2C" "In_LAMP5" "In_PVALB_Ba" "In_PVALB_Ch" "In_RELN" "In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")

for i in "${ct[@]}"; do
  Rscript ~/Documents/PostmortemBrain/workspace/scripts/RNA/08_DEanalysis/DESeq2_RELN/03_1_DESeq2-pcFromCorrectedDataOnAllCellTypes.R "$i"
done

# Error in (function (file = if (onefile) "Rplots.pdf" else "Rplot%03d.pdf",  : 
#   cannot open file 'Rplots.pdf'
# Calls: hist ... hist.default -> plot -> plot.histogram -> dev.hold -> <Anonymous>
# Execution halted