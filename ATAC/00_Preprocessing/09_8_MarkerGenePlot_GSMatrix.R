## Concatenate gene score matrices of different chromosomes
## and create marker gene plot for publication

# run with conda environment sc-atac

library(ArchR)
library(dplyr)
library(data.table)
library(viridis)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"
atacdir <- paste0(basedir, "scripts/ATAC/")

#mg_version <- "marker-v1"
#mg_version <- "marker-v2"
mg_version <- "marker-v3"

mg_paperAnna <- c("FLT1", "SYNE2", "COBLL1", # Endothelial
            "APBB1IP", "C3", "P2RY12", #Micorglia
            "MBP", "MOBP", "RNF220", #Oligodendrocytes
            "OLIG1", "OLIG2", "PDGFRA", # OPC
            "AQP4", "GFAP", "CLU", # Astrocytes
            "GAD1", "GAD2", "SLC32A1", # Inhibitory neurons
            "SATB2", "SLC17A6", "SLC17A7" # Excitatory neurons
            )
mg_paperAnna_group <- c(rep("Endothelial", 3),
                rep("Microglia",3), #Micorglia
                rep("Oligo.", 3),#Oligodendrocytes
                rep("OPC", 3), # OPC
                rep("Astrocytes", 3), # Astrocytes
                rep("Inhibitory N.", 3), # Inhibitory neurons
                rep("Excitatory N.", 3) # Excitatory neurons
                )

names(mg_paperAnna_group) <- mg_paperAnna

# 1. Read gene score matrices for the different chromosomes
p <- paste0("GeneScoreMatrix_chr.*.rds")
files <- list.files(path = paste0(atacdir, "ArchRSubset_FinalFiltering/GeneScoreMatrices/"),
            pattern = p,
            full.names = FALSE)
chrom_list <- list()
for (f in files){
    chrom_mat <- readRDS(paste0(atacdir, "ArchRSubset_FinalFiltering/GeneScoreMatrices/",f))
    index_genes <- c(mg_paperAnna, rownames(chrom_mat)[1:2]) # at least 2 genes needed so that trnsformation from sparse to regular matrix works
    print(rownames(chrom_mat)[rownames(chrom_mat) %in% mg_paperAnna])
    chrom_mat <- as.matrix(chrom_mat[rownames(chrom_mat) %in% index_genes,])
    chrom_list[[f]] <- chrom_mat
}
chrom_all <- do.call(rbind, unname(chrom_list))
# subset matrix to actual marker genes --> only 19 out of 21 found
chrom_all_mg <- chrom_all[rownames(chrom_all) %in% mg_paperAnna,]

# 2. Load archr object (e.g. to get cell type labels)
projPM <- loadArchRProject(path = "ArchRSubset_FinalFiltering/")
# relabel In_LAMP5 to In_RELN
projPM@cellColData$celltypes_refined[projPM@cellColData$celltypes_refined == "In_LAMP5"] <- "In_RELN"
#colnames(chrom_all) <- projPM@cellColData$celltypes_refined
# add column indicating broad cell type
tmp <- projPM@cellColData
tmp$broad_celltype <- tmp$celltypes_refined
tmp$broad_celltype[tmp$broad_celltype %in% c("Astro_FB", "Astro_PP")] <- "Astrocytes"
tmp$broad_celltype[startsWith(tmp$broad_celltype, "Exc")] <- "ExcitatoryNeurons"
tmp$broad_celltype[startsWith(tmp$broad_celltype, "In")] <- "InhibitoryNeurons"
projPM <- addCellColData(
  ArchRProj = projPM,
  data = tmp$broad_celltype,
  name = "broad_celltype",
  cells = rownames(tmp)
)

# 3. Get metrics for Marker Genes
#celltypes <- unique(projPM@cellColData$broad_celltype)
#celltypes_cells <- projPM@cellColData[colnames(chrom_all_mg),]$broad_celltype
celltypes <- unique(projPM@cellColData$celltypes_refined)
celltypes_cells <- projPM@cellColData[colnames(chrom_all_mg),]$celltypes_refined

get_metrics_ct <- function(ct, mat = chrom_all_mg){
    indices_col <- which(celltypes_cells == ct)
    mat_ct <- mat[,indices_col]

    frac_thresh <- apply(mat_ct, 1, function(x) sum(x > 1.0))/ncol(mat_ct)
    mean_acc <- apply(mat_ct, 1, function(x) mean(x))

    df <- data.frame(gene = names(frac_thresh),
                    frac_acc = frac_thresh,
                    mean_acc = mean_acc)
    return(df)
}

metric_list <- lapply(celltypes, get_metrics_ct)
names(metric_list) <- celltypes
metric_df <- bind_rows(metric_list, .id="celltype")

# 4. Actual Dot Plot
#metric_df$gene <- factor(metric_df$gene, levels=mg_paperAnna)
metric_df$group <- mg_paperAnna_group[metric_df$gene]
#ggplot(metric_df, aes(x=gene, y=celltype, size=frac_acc, col=log(mean_acc+1))) +
ggplot(metric_df, aes(x=gene, y=celltype, size=frac_acc, col=mean_acc)) +
    geom_point() +
    scale_x_discrete(drop=FALSE,
						guide = guide_axis(angle = 45)) +
    scale_color_viridis(name = "Mean gene score", option = "A")+
    scale_size_continuous(labels = scales::percent,
                        name = "Fraction of cells with\ngene score > 1.0") +
    facet_grid(.~group, scales = "free", switch = "x", space = "free_x") + 
    theme_light() +
	theme(
			strip.text.x = element_text(size = 12),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x = element_text(size = 12),
			axis.text.y = element_text(size = 12),
			legend.text = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
		)
# ggsave(filename = "ArchRSubset_FinalFiltering/Plots/DotPlot_MarkerGenes.pdf",
# 	width=12, height=8)
ggsave(filename = paste0("ArchRSubset_FinalFiltering/Plots/09_8_DotPlot_MarkerGenes_allCelltypes_", mg_version, ".pdf"),
	width=12, height=8)