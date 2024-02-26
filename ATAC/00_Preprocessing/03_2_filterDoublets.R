## Filter out doublets according to previously calculated doublet scores
## and filter out samples with low quality

# run with conda environment sc-atac

# 0. load libraries
library(ArchR)
library(dplyr)
library(RColorBrewer)
library(viridis)

addArchRThreads(20)

# 1. Read data

projPM1 <- loadArchRProject(path = "ArchROutput/")

# 2a. Collecting sample statistics

col_pal <- colorRampPalette(brewer.pal(name="Dark2", n = 8))(16)

s <- as.data.frame(projPM1@sampleColData[c("ID","X6.Batch")])
c <- as.data.frame(projPM1@cellColData[c("Sample","TSSEnrichment","nFrags","DoubletEnrichment")])
sc <- inner_join(c, s, by=c("Sample" = "ID"))
sc$X6.Batch <- as.factor(sc$X6.Batch)
sc$Sample <- factor(sc$Sample, levels=unique(sc[order(sc$X6.Batch),]$Sample), ordered=TRUE)


# 2b. Boxplot DoubletEnrichment scores before filtering

ggplot(sc, aes(x=Sample, y=DoubletEnrichment, fill=X6.Batch)) +
        geom_boxplot() +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-DoubletEnrichment-woFiltering-boxplot.pdf", width = 12, height = 8)

# 3. Filter out doublets
projPM2 <- filterDoublets(projPM1,
	filterRatio = 2.5)

# 4a. Collecting new sample statistics
s <- as.data.frame(projPM2@sampleColData[c("ID","X6.Batch")])
c <- as.data.frame(projPM2@cellColData[c("Sample","DoubletEnrichment")])
sc <- inner_join(c, s, by=c("Sample" = "ID"))
sc$X6.Batch <- as.factor(sc$X6.Batch)
sc$Sample <- factor(sc$Sample, levels=unique(sc[order(sc$X6.Batch),]$Sample), ordered=TRUE)

# 4b. Boxplot DoubletEnrichment scores after filtering

ggplot(sc, aes(x=Sample, y=DoubletEnrichment, fill=X6.Batch)) +
        geom_boxplot() +
        scale_fill_manual(values = col_pal) +
        scale_x_discrete(guide = guide_axis(angle = 90))
        #scale_fill_viridis(discrete=TRUE)
ggsave("ArchROutput/Plots/QC-Sample-DoubletEnrichment-wFiltering2.5ratio-boxplot.pdf", width = 12, height = 8)

## 5. Remove cells with DoubletEnrichment above 8
#
#projPM2@cellColData <- projPM2@cellColData[projPM2@cellColData$DoubletEnrichment < 5,,drop=FALSE]
#
## 6a. Collecting new sample statistics
#s <- as.data.frame(projPM2@sampleColData[c("ID","X6.Batch")])
#c <- as.data.frame(projPM2@cellColData[c("Sample","DoubletEnrichment")])
#sc <- inner_join(c, s, by=c("Sample" = "ID"))
#sc$X6.Batch <- as.factor(sc$X6.Batch)
#sc$Sample <- factor(sc$Sample, levels=unique(sc[order(sc$X6.Batch),]$Sample), ordered=TRUE)
#
## 6b. Boxplot DoubletEnrichment scores after filtering
#
#ggplot(sc, aes(x=Sample, y=DoubletEnrichment, fill=X6.Batch)) +
#        geom_boxplot() +
#        scale_fill_manual(values = col_pal) +
#        scale_x_discrete(guide = guide_axis(angle = 90))
#        #scale_fill_viridis(discrete=TRUE)
#ggsave("ArchROutput/Plots/QC-Sample-DoubletEnrichment-wFiltering5DEnrich-boxplot.pdf", width = 12, height = 8)

# 7. Filter out sample SU789

projPM2@cellColData <- projPM2@cellColData[projPM2@cellColData$Sample != "SU789_PFC_ATAC",,drop=FALSE]


# 8. Save filtered project
saveArchRProject(ArchRProj = projPM2)
