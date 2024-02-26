## Plot list of known marker genes with 
## UMAP for each gene with imputed gene scores as colour

# run with conda environment sc-atac

library(ArchR)

# 1. Load ArchR project with ATACseq data
projPM <- loadArchRProject(path = "ArchROutput/")

# 2. List of marker genes across cell types
# FeatureName (AC012405.1,TPRS1,PVRL3,TLE) does not exist!
markerGenes  <- c(
    "AQP4", "ATP1A2", "GJA1", "CLU", #Astrocytes
    "GFAP", "WDR49", "SLC38A1", "LINC01411", "TJP2", "PFKFB2", #Astrocytes fibrous
    "APOE", "ALDH1L1", "ACSBG1", "SLC1A2", "SLC39A12", "TEAD1", #Astrocytes protoplasmic
    "GLUL", "SOX9", "NDRG2", "GFAP", "ALDH1A1", "VIM", #Astrocytes Nagy
    "SYNE2","EPAS1", "ABCB1", "FLT1", "SLC7A5", "FN1", "COBLL1", "PRKCH", "VWF", "HLA-E", #Endothelial
    "CLDN5", "VTN", # Endotheliall Nagy
    "SNAP25", "STMN2", "RBFOX3", #Neurons Nagy
    "SATB2", "SLC17A7", "SLC17A6", # Ex Neurons Nagy
    "TESPA1", "CUX2", "LINC01378", "GLRA3", "LAMP5", #Ex_L2-3
    "PRSS12", #Ex_L3-5
    "THEMIS", "LINC00343", #Ex_L4-6_THEMIS_LINC00343
    "TLE4", "ETV1", "HTR2C", #Ex_L5-6
    "CLSTN2", "HS3ST4", #L4-6
    "RORB", #expressed in L4, but not in L5-6
    "TSHZ2", "CPNE4", "DCC", #L4
    "RORA", "NFIB", "UNC5C", #Distinguish In_PVALB_Ch from IN_PVALB_Ba
    "GAD2", "SLC32A1", #In Neurons Nagy
    "GAD1", "LAMP5", "GAD2", #In_LAMP5
    "PVALB", #In_PVALB
    "SST", #In_SST
    "VIP", "CALB2", #In_VIP
    "RELN", #In_RELN
    "CARTPT", #L2-3 Nagy
    "THSD7A", #L2-4 Nagy
    "RASGRF2", #L2-6 Nagy
    "RORB", #L4-5 Nagy
    "GRIK4", #L4-6 Nagy
    "KCNK2", "SULF2", "PCP4", "FEZF2", #L5 Nagy
    "TOX", "RPRM", "RXFP1", "FOXP2", #L5-6 Nagy
    "SYT6", "OPRK1", "NR4A2", "SYNPR", "NTNG2", "ADRA2A", # L6 Nagy
    "C3", "DOCK8", "APBB1IP", "P2RY12", "MEF2A", "SORL1", "PLXDC2", #Microglia Anna
    "SPI1", "MRC1", "TMEM119", "CX3CR1", # Macrophages/Migroglia Nagy
    "CTNNA3", "ST18", "RNF220", "MBP", "SLC44A1", "PLP1", #Oligodendrocytes
    "MAG", "MOG", "MOBP", "MBP", #Olidodendrocytes Nagy
    "PDGFRA", "VCAN", "MEGF11", "LHFPL3", "TNR", "PCDH15", "PTPRZ1", "DSCAM", #OPCs
    "PTGDS", "OLIG2", "OLIG1" # OPCs Nagy
  )

# 3.Plot UMAP per known marker gene of cell type using imputed weightss
p <- plotEmbedding(
    ArchRProj = projPM, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projPM)
)

# 4. Save plots as PDF
plotPDF(plotList = p, 
    name = "09_2_Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = projPM, 
    addDOC = FALSE, width = 5, height = 5)
