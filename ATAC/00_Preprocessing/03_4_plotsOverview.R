## Plotting QC statistics one more time to have a final check

# INFO: use conda environment sc-atac to run this

# 0. load libraries
library(ArchR)
library(dplyr)
library(RColorBrewer)
library(viridis)

addArchRThreads(20)

# 1. Read data

projPM1 <- loadArchRProject(path = "ArchROutput/")

# 2. plotting sample statistics

col_pal <- colorRampPalette(brewer.pal(name="Dark2", n = 8))(16)

s <- as.data.frame(projPM1@sampleColData[c("ID","X6.Batch")])
c <- as.data.frame(projPM1@cellColData[c("Sample","TSSEnrichment","nFrags","DoubletEnrichment")])
sc <- inner_join(c, s, by=c("Sample" = "ID")) 
sc$X6.Batch <- as.factor(sc$X6.Batch)
sc$Sample <- factor(sc$Sample, levels=unique(sc[order(sc$X6.Batch),]$Sample), ordered=TRUE)


# violin plot for TSS enrichment score 
ggplot(sc, aes(x=Sample, y=TSSEnrichment, fill=X6.Batch)) + 
	geom_violin() +
	scale_fill_manual(values = col_pal) +
	scale_x_discrete(guide = guide_axis(angle = 90))
	#scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-TSSEnrich-ggplot.pdf", width = 12, height = 8)


# violin plott for number of fragments
ggplot(sc, aes(x=Sample, y=log10(nFrags), fill=X6.Batch)) +
        geom_violin() +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-nFrags-ggplot.pdf", width = 12, height = 8)


# violin plott for doublet enrichment
ggplot(sc, aes(x=Sample, y=DoubletEnrichment, fill=X6.Batch)) +
        geom_violin() +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-DoubletEnrichment-ggplot.pdf", width = 12, height = 8)


# boxplot plot for TSS enrichment score 
ggplot(sc, aes(x=Sample, y=TSSEnrichment, fill=X6.Batch)) + 
	geom_boxplot() +
	scale_fill_manual(values = col_pal) +
	scale_x_discrete(guide = guide_axis(angle = 90))
	#scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-TSSEnrich-boxplot.pdf", width = 12, height = 8)

# boxplot plott for number of fragments
ggplot(sc, aes(x=Sample, y=log10(nFrags), fill=X6.Batch)) +
        geom_boxplot() +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-nFrags-boxplot.pdf", width = 12, height = 8)

# boxplot plott for Doublet Enrichment
ggplot(sc, aes(x=Sample, y=DoubletEnrichment, fill=X6.Batch)) +
        geom_boxplot() +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-DoubletEnrichment-boxplot.pdf", width = 12, height = 8)


# violin/boxplot plot for TSS enrichment score 
ggplot(sc, aes(x=Sample, y=TSSEnrichment, fill=X6.Batch)) + 
	geom_violin() +
	geom_boxplot(width=0.1, color = "grey", alpha=0.2) +
	scale_fill_manual(values = col_pal) +
	scale_x_discrete(guide = guide_axis(angle = 90))
	#scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-TSSEnrich-violinboxplot.pdf", width = 12, height = 8)

# boxplot plott for number of fragments
ggplot(sc, aes(x=Sample, y=log10(nFrags), fill=X6.Batch)) +
        geom_violin() +
	geom_boxplot(width=0.1, color="grey", alpha = 0.2) +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-nFrags-violinboxplot.pdf", width = 12, height = 8)

# boxplot plott for Doublet Enrichment
ggplot(sc, aes(x=Sample, y=DoubletEnrichment, fill=X6.Batch)) +
        geom_violin() +
	geom_boxplot(width=0.1, color="grey", alpha = 0.2) +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-DoubletEnrichment-violinboxplot.pdf", width = 12, height = 8)


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
