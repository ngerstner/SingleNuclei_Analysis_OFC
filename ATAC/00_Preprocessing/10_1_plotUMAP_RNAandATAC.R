## Plot RNA and ATAC UMAPs next to each other
## and number of nuclei per celltype for both modalities

# run with conda environment sc-atac

library(ArchR)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggpattern)

# define pathes
basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"

# paper or doctoral thesis
# mode <- "paper"
mode <- "thesis"

# read RNA data
rna <- readRDS("../../scanpy_adata/adata_labelTransfer_celltypes_samplesFilt_AnnaAnnotation_RELN_woObs.rds")
# subset RNA data to exclude Exc_7 and Exc_20
Idents(object = rna)
Idents(object = rna) <- "ctAnna_r1"
rna_filt <- subset(x = rna, idents = c("Exc_7", "Exc_20"), invert = TRUE)
rna <- rna_filt

# read ATAC data
#projPM <- loadArchRProject(path = "ArchROutput/")
# ArchR object after final filtering
projPM <- loadArchRProject(path = paste0(basedir,"scripts/ATAC/ArchRSubset_FinalFiltering/"))
# relabel In_LAMP5 to In_RELN
projPM@cellColData$celltypes_refined[projPM@cellColData$celltypes_refined == "In_LAMP5"] <- "In_RELN"

#celltypes <- unique(rna@meta.data$Anna_celltypes)
# refined cell types in rna object
celltypes <- unique(rna@meta.data$ctAnna_r1)


if (mode == "paper"){
	# colorPal <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
	#                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")
	# color palette without Exc_7 and Exc_20 (4 and 5)
	colorPal <- c("1"="#D51F26","2"="#272E6A","3"="#208A42", "6"="#FEE500", "4"="#89288F", "7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
				"10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

	names(colorPal) <- sort(celltypes)
} else if (mode == "thesis") {
   # color palette doctoral thesis
	colorPal <- readRDS(paste0(basedir, "/scripts/ColorPalette_Thesis.rds"))
}

# 1. Plot UMAPs
# plot RNA
rna_df <- as.data.frame(rna[["umap"]]@cell.embeddings)
#rna_df$celltype <- rna@meta.data$Anna_celltypes
rna_df$celltype <- rna@meta.data$ctAnna_r1

ggplot(rna_df, aes(x=umap_1, y=umap_2, col = celltype)) +
	geom_point(size = 0.5, stroke = 0, shape = 16) +
	scale_color_manual(name = "Cell type",
			breaks=names(colorPal),
			labels = names(colorPal),
			values = colorPal) +
	theme_light() +
	theme(text = element_text(size=24),
			panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
	xlab("UMAP 1") +
	ylab("UMAP 2") +
	ggtitle("snRNA-seq") +
	guides(colour = guide_legend(override.aes = list(size=10)))
if (mode == "paper"){
	# ggsave(filename = "ArchRSubset_FinalFiltering/Plots/UMAP_RNAseq_refined.png",
	# 	width = 15, height = 12)
	ggsave(filename = "ArchRSubset_FinalFiltering/Plots/UMAP_RNAseq_refined_woExc7and20.png",
		width = 15, height = 12)
} else if (mode == "thesis"){
	ggsave(filename = "ArchRSubset_FinalFiltering/Plots/UMAP_RNAseq_refined_woExc7and20_thesis.png",
		width = 15, height = 12)
}


# plot ATAC
atac_df <- projPM@embeddings$UMAP$df
colnames(atac_df) <- c("umap1", "umap2")
atac_df$celltype <- as.vector(projPM@cellColData$celltypes_refined)

ggplot(atac_df, aes(x=umap1, y = umap2, col = celltype)) +
	geom_point(size = 0.5, stroke = 0, shape = 16) +
	scale_color_manual(name = "Cell type",
			breaks=names(colorPal),
			labels = names(colorPal),
			values = colorPal) +
	theme_light() +
	theme(text = element_text(size = 24),
			panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
	xlab("UMAP 1") +
	ylab("UMAP 2") +
	ggtitle("snATAC-seq") +
	guides(colour = guide_legend(override.aes = list(size=10)))
if (mode == "paper"){
	ggsave(filename = "ArchRSubset_FinalFiltering/Plots/UMAP_ATACseq_finalFiltering.png",
	width=15, height=12)
} else if (mode == "thesis"){
	ggsave(filename = "ArchRSubset_FinalFiltering/Plots/UMAP_ATACseq_finalFiltering_thesis.png",
	width=15, height=12)
}


# 2. Plot number of nuclei per celltype for both data modalities

n_rna <- table(rna@meta.data$ctAnna_r1)
n_atac <- table(projPM@cellColData$celltypes_refined)
n_df <- as.data.frame(n_rna) %>%
	full_join(as.data.frame(n_atac), by="Var1") %>%
	dplyr::rename(RNA = Freq.x, ATAC = Freq.y, celltype = Var1) 
write.csv(n_df, paste0(basedir, "tables/NumberOfNuclei_perCelltype.csv"))
n_df <- n_df %>%
	tidyr::pivot_longer(cols=RNA:ATAC) %>%
	mutate(name = factor(name, levels = c("RNA", "ATAC")))

ggplot(data = n_df, aes(x = celltype, y=value, fill = celltype, pattern = name)) +
  	geom_bar_pattern(stat="identity",
					position = position_dodge(preserve = "single"),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(name = "Cell type", 
  			breaks=names(colorPal),
			labels = names(colorPal),
			values = colorPal) +
  scale_pattern_manual(values = c(ATAC = "stripe", RNA = "none")) +
  scale_x_discrete(drop=FALSE,
						guide = guide_axis(angle = 45)) +
  labs(x = "Cell type", y = "Number of Nuclei", pattern = "Data modality") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
	theme_light() +
	theme(
			strip.text.x = element_text(size = 15),
			axis.title.x = element_blank(),
			axis.title.y = element_text(size = 15),
			axis.text.x = element_text(size = 12),
			axis.text.y = element_text(size = 12),
			legend.text = element_text(size = 12)
		)
ggsave(filename = "ArchRSubset_FinalFiltering/Plots/NumberOfNuclei_PerCelltype_RNAandATAC_finalFiltering.pdf",
	width=12, height=8)


# 3. Statistics 

# 3.1 Number of nuclei per donor
# RNA
tmp_RNA <- table(rna@meta.data$sample)
df_RNA <- data.frame("individual" = names(tmp_RNA),
	"count"=tmp_RNA)
df_RNA$individual <- str_replace(df_RNA$individual, "_PFC_RNA", "")
mean(df_RNA$count.Freq)
# [1] 9046.506
min(df_RNA$count.Freq)
# [1] 3895
max(df_RNA$count.Freq)
# [1] 15693

# ATAC
tmp <- table(projPM@cellColData$Sample)
df_ATAC <- data.frame("individual" = names(tmp),
	"count"=tmp) 
df_ATAC$individual <- str_replace(df_ATAC$individual, "_PFC_ATAC", "")
mean(df_ATAC$count.Freq)
#[1] 4438.211
min(df_ATAC$count.Freq)
#[1] 982
max(df_ATAC$count.Freq)
#[1] 8707

df_comb <- full_join(df_RNA, df_ATAC, by="individual") %>%
	rename(RNA = count.Freq.x, ATAC=count.Freq.y) %>%
	select(individual, RNA, ATAC) 
write.csv(df_comb, paste0(basedir, "tables/NumberOfNuclei_perIndividual.csv"))

# 3.2 Number of UMI counts/fragments per nuclei
# RNA
mean(rna@meta.data$n_counts)
#[1] 4884.495
median(rna@meta.data$n_counts)
#[1] 3887
min(rna@meta.data$n_counts)
#[1] 500
max(rna@meta.data$n_counts)
#[1] 14208

# ATAC
mean(projPM@cellColData$nFrags)
#[1] 8991.813
median(projPM@cellColData$nFrags)
#[1] 7071
min(projPM@cellColData$nFrags)
#[1] 1000
max(projPM@cellColData$nFrags)
#[1] 96319

# 3.3 Number of genes per nucleus

#RNA 
# RNA
mean(rna@meta.data$n_genes)
#[1] 2339.062
median(rna@meta.data$n_genes)
#[1] 2205
min(rna@meta.data$n_genes)
#[1] 301
max(rna@meta.data$n_genes)
#[1] 5817