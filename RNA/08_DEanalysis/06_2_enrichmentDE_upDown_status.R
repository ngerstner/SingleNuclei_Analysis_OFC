## Enrichment analysis (Pathway, GO terms, Diseases) for 250 most
## up- and downregulated genes per celltype for PRS extreme groups

## 0. Load Packages

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(pheatmap)
library(RColorBrewer)
library(treemap)
library(gridpattern)
library(grid)


## 1. Setup

# define pathes to data and folders
base_workspace <- '~/Documents/PostmortemBrain/workspace/'
basedir <- paste0(base_workspace, 'scripts/RNA/08_DEanalysis/DESeq2_RELN/')

# type of test: Likelihood Ratio Test or Wald test
test <- "Wald"
# test <- "LRT" 


## 2. Read in DE genes
# read all genes to have them as background lists and use first n genes
p <- paste0("X6.Batch_", test, "_pcFromCorrectedDataAllCellTypes.csv")
files <- list.files(path = paste0(basedir, "tables/DESeq2_", test, "/"),
                    pattern = p,
                    full.names = FALSE)
files <- files[!str_detect(files, "metadata")]
# read in the diff genes of different cell types
de_up <- list()
de_down <- list()
for (f in files){
  de <- read.csv(paste0(basedir, "tables/DESeq2_", test, "/", f))
  ct <- sub("(\\w*)_~.*","\\1",f)
  tmp_up <- de[de$log2FoldChange > 0, ]
  tmp_down <- de[de$log2FoldChange < 0, ]
  de_up[[ct]] <- tmp_up[order(tmp_up$padj),]
  de_down[[ct]] <- tmp_down[order(tmp_down$padj),]
}
celltypes <- names(de_up)


## 3. Enrichment analyses
enr <- list()
for (d in c("up", "down")){
  de_list <- if (d == "up") de_up else de_down
  enr[[d]] <- list()
  enr[[d]][["go_enrich"]] <- list()
  enr[[d]][["disease_enrich"]] <- list()
  enr[[d]][["kegg_enrich"]] <- list()
  
  for (ct in names(de_list)){
    # entrez_genes <-mapIds(org.Hs.eg.db, keys = de_ct[[ct]]$gene,
    entrez_genes <-mapIds(org.Hs.eg.db, keys = de_list[[ct]][1:250,]$gene,
                          keytype="ENSEMBL", column = "ENTREZID")
    entrez_bg <- mapIds(org.Hs.eg.db, keys = de_list[[ct]]$gene,
                        keytype="ENSEMBL", column = "ENTREZID")
    
    # GO enrichment
    ego <- enrichGO(gene          = entrez_genes,
                    universe      = entrez_bg,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    ego_df <- ego@result
    enr[[d]][["go_enrich"]][[ct]] <- ego_df
    
    ego_df[1:25,] %>%
    mutate(Description = fct_reorder(Description, -p.adjust)) %>%
      ggplot(aes(x=Description, y=-log10(p.adjust), color=Count)) +
      geom_point(size = 4) +
      ylim(0, max(-log10(ego_df$p.adjust[1:25]))) +
      ggtitle(paste0("GO enrichment in ", ct)) +
      coord_flip()
    ggsave(filename = paste0(basedir, "plots/06_2_DEenrichments_upDown/",
                             "status/06_2_",ct,"_GOenrichment_",d,"_DE_status.pdf"))
    
    # Disease enrichment
    x <- enrichDO(gene          = entrez_genes,
                  ont           = "DO",
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  universe      = entrez_bg,
                  minGSSize     = 5,
                  maxGSSize     = 500,
                  qvalueCutoff  = 0.05,
                  readable      = FALSE)@result
    enr[[d]][["disease_enrich"]][[ct]] <- x
    x[1:20,] %>%
      mutate(Description = fct_reorder(Description, -p.adjust)) %>%
      ggplot(aes(x=Description, y=-log10(p.adjust), color=Count)) +
      geom_point(size = 4) +
      ylim(0, max(-log10(x$p.adjust[1:20]))) +
      ggtitle(paste0("Disease enrichment in ", ct)) +
      coord_flip()
    ggsave(filename = paste0(basedir, "plots/06_2_DEenrichments_upDown/", 
                             "status/06_2_",ct,"_DiseaseEnrichment_",d,"_DE_status.pdf"))
  
    # KEGG pathway enrichment
    kk <- enrichKEGG(gene         = entrez_genes,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)@result
    enr[[d]][["kegg_enrich"]][[ct]] <- kk
    kk[1:20,] %>%
      mutate(Description = fct_reorder(Description, -p.adjust)) %>%
      ggplot(aes(x=Description, y=-log10(p.adjust), color=Count)) +
      geom_point(size = 4) +
      ylim(0, max(-log10(kk$p.adjust[1:20]))) +
      ggtitle(paste0("KEGG pathway enrichment in ", ct)) +
      coord_flip()
    ggsave(filename = paste0(basedir, "plots/06_2_DEenrichments_upDown/",  
                             "status/06_2_",ct,"_KEGGEnrichment_",d,"_DE_status.pdf"))
    
    
  }
}


## 4. Combine enrichments over celltypes
df_list <- list()

for (d in c("up", "down")){
  
  df_list[[d]][["go_df"]] <- bind_rows(enr[[d]][["go_enrich"]], .id="celltype") %>%
    dplyr::select(ID, Description, p.adjust, celltype) %>%
    pivot_wider(names_from = celltype,
                values_from = p.adjust) %>%
    filter(if_any(where(is.numeric), ~ . < 0.05))
  write.csv(df_list[[d]][["go_df"]], paste0(basedir, "tables/KEGGenrichment/status/",
                                            "Regulation_", d, "_GOenrichment.csv"))
  df_list[[d]][["go_df"]] <- df_list[[d]][["go_df"]] %>% 
    mutate(across(where(is.numeric), function(x) -log10(x))) %>%
    replace(is.na(.), 0)
  
  df_list[[d]][["kegg_df"]] <- bind_rows(enr[[d]][["kegg_enrich"]], .id="celltype") %>%
    dplyr::select(ID, Description, p.adjust, celltype) %>%
    pivot_wider(names_from = celltype,
                values_from = p.adjust) %>%
    filter(if_any(where(is.numeric), ~ . < 0.05)) 
  write.csv(df_list[[d]][["kegg_df"]], paste0(basedir, "tables/KEGGenrichment/status/",
                                              "Regulation_", d, "_KEGGenrichment.csv"))
  df_list[[d]][["kegg_df"]] <- df_list[[d]][["kegg_df"]] %>%  
    mutate(across(where(is.numeric), function(x) -log10(x))) %>%
    replace(is.na(.), 0)
  
  df_list[[d]][["disease_df"]] <- bind_rows(enr[[d]][["disease_enrich"]], .id="celltype") %>%
    dplyr::select(ID, Description, p.adjust, celltype) %>%
    pivot_wider(names_from = celltype,
                values_from = p.adjust) %>%
    filter(if_any(where(is.numeric), ~ . < 0.05)) 
  write.csv(df_list[[d]][["disease_df"]], paste0(basedir, "tables/KEGGenrichment/status/",
                                                 "Regulation_", d, "_DiseaseEnrichment.csv"))
  df_list[[d]][["disease_df"]] <- df_list[[d]][["disease_df"]] %>% 
    mutate(across(where(is.numeric), function(x) -log10(x))) %>%
    replace(is.na(.), 0)
}


## 5. Plot heatmap

# plot_heatmap <- function(x, d, enr_type){
#   mat <- as.matrix(x[,3:21])
#   rownames(mat) <- pull(x, Description)
#   
#   if (enr_type == "KEGGenrichment"){
#     kegg_hierarchy <- read.csv(paste0(basedir, 'tables/KEGGenrichment/KEGGpathways_hierarchy_reformatted.csv'))
#     
#     # get hierarchical color scheme for annotations
#     kegg_hierarchy <- kegg_hierarchy %>%
#       filter(Pathway %in% rownames(mat))
#     tp <- treepalette(kegg_hierarchy, 
#                 index = c("FamilyName", "GroupName"),
#                 method = "HCL")
#     kegg_hierarchy <- kegg_hierarchy %>%
#       left_join(tp, by=c("FamilyName", "GroupName"))
#     
#     annot_colour <- list(
#       Pathway = setNames(kegg_hierarchy$HCL.color, kegg_hierarchy$Pathway)
#     )
#     pathway_col <- data.frame(Pathway = factor(kegg_hierarchy$Pathway,
#                               levels = unique(kegg_hierarchy$Pathway)))
#     row.names(pathway_col) <- kegg_hierarchy$Pathway
#     
#     # reorder matrix
#     mat <- mat[kegg_hierarchy$Pathway,]
#     mat_sig <- matrix(nrow = nrow(mat), ncol = ncol(mat))
#     mat_sig[mat>=(-log10(0.05))] <- "*"
#     mat_sig[mat<(-log10(0.05))] <- ""
#     
#     pdf(file = paste0(basedir, "plots/06_2_DEenrichments_upDown/", 
#                       "status/06_2_AllCellTypes_",enr_type,"_",d,"_DE_status.pdf"),
#         width = 12, height = 10)
#     pheatmap(as.matrix(mat),
#              display_numbers = mat_sig,
#              cluster_rows = FALSE,
#              cluster_cols = TRUE,
#              color = colorRampPalette(c("white", "#D45E60"))(100),
#              annotation_colors = annot_colour,
#              annotation_row = pathway_col
#     )
#     dev.off()
#   }
#   
# }

plot_heatmap <- function(x, d, enr_type){
  mat <- as.matrix(x[,3:21])
  rownames(mat) <- pull(x, Description)
  mat_sig <- matrix(nrow = nrow(mat), ncol = ncol(mat))
  mat_sig[mat>=(-log10(0.05))] <- "*"
  mat_sig[mat<(-log10(0.05))] <- ""
  
  pdf(file = paste0(basedir, "plots/06_2_DEenrichments_upDown/",
                    "status/06_2_AllCellTypes_",enr_type,"_",d,"_DE_status.pdf"),
      width = 12, height = 10)
  pheatmap::pheatmap(as.matrix(mat),
           display_numbers = mat_sig,
           cluster_rows = FALSE,
           cluster_cols = FALSE
  )
  dev.off()
}


for (d in c("up", "down")){
  for(et in c("GOenrichment", "DiseaseEnrichment")){
    
    if (et == "GOenrichment"){
      plot_heatmap(df_list[[d]][["go_df"]], d, et)
    } else {
      plot_heatmap(df_list[[d]][["disease_df"]], d, et)
    }
    
  }
}


plot_KEGG_heatmap <- function(x_up, x_down, enr_type){
  
  # concatenate df with up and downregulated pathways
  df_both <- bind_rows(list("up" = x_up,
                            "down" = x_down),
                       .id = "direction")
  
  # upregulated KEGG enrichment to matrix
  mat_up <- as.matrix(x_up[,3:21])
  rownames(mat_up) <- pull(x_up, Description)
  
  # downregulated KEGG enrichment to matrix
  mat_down <- as.matrix(x_down[,3:21])
  rownames(mat_down) <- pull(x_down, Description)
  
  
  # Generate color scheme representing hierarchy of KEGG pathways
  kegg_hierarchy <- read.csv(paste0(basedir, 'tables/KEGGenrichment/KEGGpathways_hierarchy_reformatted.csv'))
  
  # subset KEGG pathways to those up or downregulated
  kegg_hierarchy <- kegg_hierarchy %>%
    filter(Pathway %in% df_both$Description)
  # get colour palette with one colour per group and similar colours within family
  tp <- treepalette(kegg_hierarchy, 
                    index = c("FamilyName", "GroupName"),
                    method = "HCL")
  # add colour palette to hierarchy df
  kegg_hierarchy <- kegg_hierarchy %>%
    left_join(tp, by=c("FamilyName", "GroupName"))
  
  # get "average" colour for each family
  kegg_hierarchy <- kegg_hierarchy %>% 
    group_by(FamilyName) %>% 
    mutate(avg_color = mean_col(HCL.color))
    
  
  # get KEGG hierarchy for up and downregulated terms separately
  hierarchy_list <- list()
  hierarchy_list[["up"]] <- kegg_hierarchy %>%
    filter(Pathway %in% rownames(mat_up))
  hierarchy_list[["down"]] <- kegg_hierarchy %>%
    filter(Pathway %in% rownames(mat_down))
  hierarchy_list[["mix"]] <- kegg_hierarchy %>%
    filter(Pathway %in% rownames(mat_up) | Pathway %in% rownames(mat_down))
  
  for (d in c("up", "down")){
    # define colour palette for heatmap
    # annot_colour <- list(
    #   Group = setNames(unique(hierarchy_list[[d]]$HCL.color), unique(hierarchy_list[[d]]$GroupName)),
    #   Family = setNames(unique(hierarchy_list[[d]]$avg_color), unique(hierarchy_list[[d]]$FamilyName))
    # )
    annot_colour <- list(
        Group = setNames(unique(hierarchy_list[["mix"]]$HCL.color), unique(hierarchy_list[["mix"]]$GroupName)),
        Family = setNames(unique(hierarchy_list[["mix"]]$avg_color), unique(hierarchy_list[["mix"]]$FamilyName))
      )
    # annotation dataframe for group and family categories
    group_col <- data.frame(Group = factor(hierarchy_list[[d]]$GroupName,
                                               levels = unique(hierarchy_list[[d]]$GroupName)),
                            Family = factor(hierarchy_list[[d]]$FamilyName,
                                            levels = unique(hierarchy_list[[d]]$FamilyName)))
    row.names(group_col) <- hierarchy_list[[d]]$Pathway
      
    # reorder matrix
    mat <- if(d == "up") mat_up else mat_down
    mat <- mat[hierarchy_list[[d]]$Pathway,]
    mat_sig <- matrix(nrow = nrow(mat), ncol = ncol(mat))
    mat_sig[mat>=(-log10(0.05))] <- "*"
    mat_sig[mat<(-log10(0.05))] <- ""
    
    pdf(file = paste0(basedir, "plots/06_2_DEenrichments_upDown/", 
                      "status/06_2_AllCellTypes_",enr_type,"_",d,"_DE_status.pdf"),
        width = 12, height = 10)
    pheatmap::pheatmap(as.matrix(mat),
             display_numbers = mat_sig,
             cluster_rows = FALSE,
             cluster_cols = TRUE,
             color = colorRampPalette(c("white", "#D45E60"), bias = )(100),
             breaks = c(seq(0,5,length.out=50),
                        seq(5.2,max(max(mat_up), max(mat_down)),length.out=51)),
             annotation_colors = annot_colour,
             annotation_row = group_col,
             # fontsize = 12,
             angle_col = 45
    )
    dev.off()
  }
  
}
plot_KEGG_heatmap(df_list[["up"]][["kegg_df"]], 
                  df_list[["down"]][["kegg_df"]], 
                  "KEGG_enrichment")
# x_up <- df_list[["up"]][["kegg_df"]]
# x_down <- df_list[["down"]][["kegg_df"]]
# enr_type <- "KEGG_enrichment"


# # 6. Simplify redundance in GO terms
# library(rrvgo)
# simMatrix <- calculateSimMatrix(df_list$up$go_df$ID,
#                                 orgdb="org.Hs.eg.db",
#                                 ont="BP",
#                                 method="Rel")
# scores <- setNames(-log10(df_list$up$go_df$Astro_FB), df_list$up$go_df$ID)
# reducedTerms <- reduceSimMatrix(simMatrix,
#                                 scores,
#                                 threshold=0.7,
#                                 orgdb="org.Hs.eg.db")
# treemapPlot(reducedTerms)
