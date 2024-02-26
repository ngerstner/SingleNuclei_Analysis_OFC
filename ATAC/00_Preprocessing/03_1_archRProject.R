## Create ArchR project and plot some statistics for QC,
## e.g. TSS enrichment, number of fragments and fragment size

# INFO: use conda environment sc-atac to run this

# 0. load libraries
library(ArchR)
library(dplyr)

# 1. general setting
addArchRThreads(threads = 20)
addArchRGenome("hg38")


# 2. generate list with all arrow files
arrowfolder <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/scripts/ATAC/ArrowFiles"

ArrowFiles <- list.files(path = arrowfolder, pattern = "\\.arrow", full.names = TRUE)


# 3. create ArchRProject

projPM1 <- ArchRProject(
	ArrowFiles = ArrowFiles,
	outputDirectory = "ArchROutput",
	copyArrows = TRUE
)

print(projPM1)
print(paste0("Memory Size = ", round(object.size(projPM1) / 10^6, 3), " MB"))
print(getAvailableMatrices(projPM1))


# 4. add phenotype data

pheno <- read.csv("/psycl/g/mpsngs/HiSeq_Helmholtz/20210324_Anna_Froehlich_10X_RNAseq/phenotype_clean.csv", header = TRUE)
pheno <- pheno %>%
	mutate(ID = paste0("SU",SU.Number,"_PFC_ATAC"))

projPM1@sampleColData$ID <- rownames(projPM1@sampleColData)

projPM1@sampleColData <- merge(projPM1@sampleColData, pheno, by = "ID")
rownames(projPM1@sampleColData) <- projPM1@sampleColData$ID

# 5. plotting sample statistics

# ridge plot for each sample for TSS enrichment scores
p1 <- plotGroups(
	ArchRProj = projPM1,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "TSSEnrichment",
	plotAs = "ridges"
)

# violin plot for each sample for TSS enrichment scores
p2 <- plotGroups(
	ArchRProj = projPM1,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "TSSEnrichment",
	plotAs = "violin",
	alpha = 0.4,
	addBoxPlot = TRUE
)

# ridge plot for each sample for log10(unique nuclear fragments)
p3 <- plotGroups(
	ArchRProj = projPM1,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "log10(nFrags)",
	plotAs = "ridges"
)

# violin plot for each sample for log10(unique nuclear fragments)
p4 <- plotGroups(
	ArchRProj = projPM1,
	groupBy = "Sample",
	colorBy = "cellColData",
	name = "log10(nFrags)",
	plotAs = "violin",
	alpha = 0.4,
	addBoxPlot = TRUE
)

# save plots in PDF file
plotPDF(p1,p2,p3,p4, 
	name = "QC-Sample-Statistics.pdf",
	ArchRProj = projPM1,
	addDOC = FALSE,
	width = 12, 
	height = 12
)


# 6. plotting sample fragment size distribution and TSS enrichment profiles

# fragment size distribution
pf1 <- plotFragmentSizes(
	ArchRProj = projPM1
)

plotPDF(pf1,
	name = "QC-Sample-FragSizes.pdf",
	ArchRProject = projPM1,
	addDOC = FALSE,
	width = 12,
	height = 12
)

# TSS enrichment profiles
pf2 <- plotTSSEnrichment(
	ArchRProj = projPM1
)

# save plots in PDF file
plotPDF(pf1, pf2,
	name = "QC-Sample-FragSizes-TSSProfile.pdf",
	ArchRProj = projPM1,
	addDOC = FALSE, 
	width = 12,
	height = 12
)


# 7. save ArchRProject

saveArchRProject(
	ArchRProj = projPM1
)
