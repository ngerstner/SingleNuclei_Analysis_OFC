## Plot distributions of number of nuclei of individuals, batches etc. per cluster
## e.g. heatmap with donor vs. cell type --> color indicates number of nuclei

# run with conda environment DESeq2

## 0. Load packages
library(dplyr)
library(ggplot2)
library(tibble)
library(ComplexHeatmap)

basedir <- "/psycl/g/mpsagbinder/mgp/workspace/SingleNuc_PostmortemBrain/"

# 1. Read RNA and ATAC cell type proportions

ct_RNA <- read.csv(paste0(basedir, "tables/distrCelltypePerSample_RELN.csv")) %>%
        remove_rownames %>% 
        column_to_rownames(var="sample")
colnames(ct_RNA) <- gsub("\\.", "-", colnames(ct_RNA))
rownames(ct_RNA) <- gsub("_RNA", "", rownames(ct_RNA))
ct_RNA <- as.matrix(ct_RNA)
head(ct_RNA)

ct_ATAC <- read.csv(paste0(basedir, "tables/distrCelltypePerSample_ATAC_RELN.csv")) %>%
        remove_rownames %>%
        column_to_rownames(var="X")
ct_ATAC <- t(as.matrix(ct_ATAC))
rownames(ct_ATAC) <- gsub("_ATAC", "", rownames(ct_ATAC))
head(ct_ATAC)

# 2. Read phenotype data to annotate status in heatmap

pheno <- read.csv(paste0(basedir, "phenotype/phenotype_clean.csv"))
pheno$ID <- paste0("SU", pheno$SU.Number, "_PFC")
rownames(pheno) <- pheno$ID

id_case <- pheno$ID[pheno$Status == 1]
id_control <- pheno$ID[pheno$Status == 0]

# 3. Plot heatmap for RNA and ATAC distributions

pheno_RNA <- pheno[rownames(ct_RNA),]
column_ha <- HeatmapAnnotation(Status = as.factor(pheno_RNA$Status),
            col = list(Status = c("0" = "#EFEE68", "1" = "#48A4B6")))
pdf(paste0(basedir, "scripts/ATAC/ArchRSubset_FinalFiltering/Plots/10_6_CellTypeProportions/10_6_RNA_DistrCelltypePerSample.pdf"),
    width=15, height=8)
Heatmap(t(ct_RNA), cluster_rows = FALSE, top_annotation = column_ha, #, cluster_columns = FALSE)
        col = RColorBrewer::brewer.pal(name = "Oranges", n = 9))
dev.off()

pheno_ATAC <- pheno[rownames(ct_ATAC),]
column_ha <- HeatmapAnnotation(Status = as.factor(pheno_ATAC$Status),
            col = list(Status = c("0" = "#EFEE68", "1" = "#48A4B6")))
pdf(paste0(basedir, "scripts/ATAC/ArchRSubset_FinalFiltering/Plots/10_6_CellTypeProportions/10_6_ATAC_DistrCelltypePerSample.pdf"),
    width=15, height=8)
Heatmap(t(ct_ATAC), cluster_rows = FALSE, top_annotation = column_ha, #, cluster_columns = FALSE)
        col = RColorBrewer::brewer.pal(name = "Oranges", n = 9))
dev.off()


# 4. Correlation between proportions in RNA and ATAC

intersect_ct <- intersect(colnames(ct_RNA), colnames(ct_ATAC))
intersect_ct

intersect_sample <- intersect(rownames(ct_RNA), rownames(ct_ATAC))
intersect_sample

ct_RNA_common <- ct_RNA[intersect_sample,intersect_ct]
ct_ATAC_common <- ct_ATAC[intersect_sample,intersect_ct]

cor_mod <- cor(t(ct_RNA_common), t(ct_ATAC_common))

pdf(paste0(basedir, "scripts/ATAC/ArchRSubset_FinalFiltering/Plots/10_6_CellTypeProportions/10_6_Correlation_CelltypePorportions_RNAandATAC.pdf"),
    width=15, height=15)
Heatmap(cor_mod, cluster_rows = FALSE, cluster_columns = FALSE,
        col = RColorBrewer::brewer.pal(name = "Blues", n = 9))
dev.off()


cor_diag <- diag(cor_mod)
pdf(paste0(basedir, "scripts/ATAC/ArchRSubset_FinalFiltering/Plots/10_6_CellTypeProportions/10_6_Correlation_CelltypePorportions_RNAandATAC_HistogramDiagonal.pdf"),
    width=10, height=8)
hist(cor_diag, 
main="Histogram of correlations between celltype proportions in RNAseq and ATACseq data",
xlab="Correlation of celltype proportions",
xlim=c(0,1),
col="darkblue")
dev.off()

cor_diag <- as.data.frame(cor_diag) %>%
                rownames_to_column(var = "sample")
ggplot(cor_diag, aes(x=cor_diag)) + 
  geom_histogram(color="black", fill="darkblue") +
  geom_vline(xintercept = median(cor_diag$cor_diag), color="darkred", linetype="dashed") +
  xlab("Correlation of celltype proportions") +
  ylab("Frequency") +
  xlim(0,1) +
  theme_light() +
	theme(
			strip.text.x = element_text(size = 15),
			axis.title.x = element_text(size = 15),
			axis.title.y = element_text(size = 15),
			axis.text.x = element_text(size = 12),
			axis.text.y = element_text(size = 12),
			legend.text = element_text(size = 12)
		)
ggsave(paste0(basedir, "scripts/ATAC/ArchRSubset_FinalFiltering/Plots/10_6_CellTypeProportions/10_6_Correlation_CelltypePorportions_RNAandATAC_HistogramDiagonal_ggplot.pdf"),
    width=10, height=7)

# 5. Two proportions z-test to test for differences in distribution

plot_histogram_prop <- function(vec1, vec2, name1, name2, res_pval){

        plot_df <- data.frame("Group" = c(rep(name1, length(vec1)), rep(name2, length(vec2))),
                        "Value" = c(vec1, vec2))
        # dataframe for mean value per group
        mu <- plot_df %>%
                group_by(Group) %>%
                dplyr::summarize(group_mean = mean(Value, na.rm=TRUE))
        ggplot(plot_df, aes(x=Value, fill=Group, color=Group)) +
                geom_histogram(position="identity", alpha=0.5) +
	        geom_vline(data=mu, aes(xintercept=group_mean, color=Group),
                        linetype="dashed") +
                ggtitle(paste(ct, "- t-test pvalue: ", res_pval)) +
                xlab("Cell type proportions") +
                xlim(c(0, 1)) +
                theme_light() +
                theme(
                        axis.title.x = element_text(size=15),
                        axis.title.y = element_text(size=15),
                        axis.text.x = element_text(size=12),
                        axis.text.y = element_text(size=12),
                        legend.text = element_text(size=12)
                )
                ggsave(filename = paste0(basedir, "scripts/ATAC/ArchRSubset_FinalFiltering/Plots/",
                                        "10_6_CellTypeProportions/10_6_CellTypeProportions_", ct, "_", name1, "-", name2, ".pdf"),
	                width= 8,
	                height=6)

}


pval_list_RNA <- list()
pval_list_ATAC <- list()
pval_list_modalities <- list()
pval_list_all <- list()
for (ct in intersect_ct){

        # subset IDs to those in the matrices
        id_case_RNA <- id_case[id_case %in% rownames(ct_RNA)]
        id_control_RNA <- id_control[id_control %in% rownames(ct_RNA)]
        id_case_ATAC <- id_case[id_case %in% rownames(ct_ATAC)]
        id_control_ATAC <- id_control[id_control %in% rownames(ct_ATAC)]
        id_case_modalities <- id_case[id_case %in% intersect(rownames(ct_RNA), rownames(ct_ATAC))]

        # mean(ct_RNA[id_case_RNA,ct])
        # mean(ct_ATAC[id_control_RNA,ct])
        # median(ct_RNA[id_case_RNA,ct])
        # median(ct_ATAC[id_control_RNA,ct])
        # var(ct_RNA[id_case_RNA,ct])
        # var(ct_ATAC[id_control_RNA,ct])
        res_RNA <- t.test(ct_RNA[id_case_RNA,ct], ct_RNA[id_control_RNA,ct], var.equal = TRUE)
        pval_list_RNA[[ct]] <- res_RNA$p.value
        res_ATAC <- t.test(ct_ATAC[id_case_ATAC,ct], ct_ATAC[id_control_ATAC,ct], var.equal = TRUE)
        pval_list_ATAC[[ct]] <- res_ATAC$p.value
        res_modalities <- t.test(ct_RNA[id_case_modalities,ct], ct_ATAC[id_case_modalities,ct], paired = TRUE)
        pval_list_modalities[[ct]] <- res_modalities$p.value
        res_all <- t.test(ct_RNA[intersect_sample,ct], ct_ATAC[intersect_sample,ct], paired = TRUE)
        pval_list_all[[ct]] <- res_all$p.value

        plot_histogram_prop(ct_RNA[id_case_RNA,ct], ct_RNA[id_control_RNA,ct], "RNAcases", "RNAcontrols", res_RNA$p.value)
        plot_histogram_prop(ct_ATAC[id_case_ATAC,ct], ct_ATAC[id_control_ATAC,ct], "ATACcases", "ATACcontrols", res_ATAC$p.value)
        plot_histogram_prop(ct_RNA[id_case_modalities,ct], ct_ATAC[id_case_modalities,ct], "RNAcases", "ATACcases", res_RNA$p.value)
        plot_histogram_prop(ct_RNA[intersect_sample,ct], ct_ATAC[intersect_sample,ct], "RNA", "ATAC", res_all$p.value)
}

# reshape pvalue lists to dataframes
plot_list <- function(l, name1, name2){
        df <- as.data.frame(unlist(l))
        colnames(df) <- "pvalue"
        df$celltype <- rownames(df)
        df$FDR <- p.adjust(df$pvalue, method="fdr")
        write.csv(df, paste0(basedir, "tables/TestCelltypeProportions_", name1, "-", name2, ".csv"))

        # color palette without Exc_7 and Exc_20 (4 and 5)
        colorPal <- c("1"="#D51F26","2"="#272E6A","3"="#208A42", "6"="#FEE500", "4"="#89288F", "7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

        names(colorPal) <- sort(colnames(ct_RNA))

        ggplot(df, aes(x=celltype, y=-log10(FDR), fill = celltype)) +
                geom_bar(stat = "identity", color = "black") + 
                scale_x_discrete(drop=FALSE,
                        guide = guide_axis(angle=45)) +
                scale_fill_manual(guide="none",
                        breaks = sort(df$celltype),
                        values = colorPal[sort(df$celltype)]) +
                xlab("Celltype") +
                ylab("-log10(FDR)") +
                #ggtitle(paste0("Comparison celltype proportions between ", name1, " and ", name2)) +
                theme_light() +
                theme(
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size=15),
                        axis.text.x = element_text(size=12),
                        axis.text.y = element_text(size=12),
                        legend.text = element_text(size=12)
                )
        ggsave(filename = paste0(basedir, "scripts/ATAC/ArchRSubset_FinalFiltering/Plots/",
                                        "10_6_CellTypeProportions/10_6_CellTypeProportionsPvalues_", name1, "-", name2, ".pdf"),
	        width= 10,
	        height=7)

}
plot_list(pval_list_RNA, "RNAcases", "RNAcontrols")
plot_list(pval_list_ATAC, "ATACcases", "ATACcontrols")
plot_list(pval_list_modalities, "RNAcases", "ATACcases")
plot_list(pval_list_all, "RNA", "ATAC")