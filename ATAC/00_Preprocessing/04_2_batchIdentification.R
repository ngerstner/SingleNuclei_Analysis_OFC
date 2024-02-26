## Identify batch effects by correlating covariates
## to first 10 LSI components

# run with conda environment sc-atac

library(ArchR)
library(dplyr)

addArchRThreads(20)

# 1. Read data

projPM <- loadArchRProject(path = "ArchROutput/")

# 2. Prepare data

lsi <- as.data.frame(getReducedDims(projPM))[,1:10]

cellMeta <- as.data.frame(projPM@cellColData)
sampleMeta <- as.data.frame(projPM@sampleColData)

cov <- left_join(cellMeta, sampleMeta, by=c("Sample"="ID"))


# 3. Estimate correlation/dependency between covariates and first 10 LSI components

batch_df <- data.frame("LSIcomponent" = character(), 
			"Covariate" = character(),
			"correlation" = numeric(),
			"p-value" = numeric(),
			"class" = character()) #categorical or numeric covariate
cov_test_categorical <- c("Sex", "Status", "Classification", "Hemisphere", "Mode.of.death", "Antipsychotics",
				"Smoking.status", "Main.Batch", "X6.Batch.y", "Position", "Operator")
cov_test_numerical <- c("Age", "PMI", "Brain.pH", "RIN", "Freezer.storage.time..days.", "TSSEnrichment", "nFrags", "mito_perc") 

# test categorical variables
for (l in colnames(lsi)){
	for (c in cov_test_categorical){
		anova <- aov(lsi[[l]] ~ as.factor(cov[[c]]))
		p <- summary(anova)[[1]][["Pr(>F)"]][1]
		coef <- max(abs(anova$coefficients))			
		batch_df <- rbind(batch_df, list(l, c, coef, p, "categorical"))
	}
}

# test continuous variables
for (l in colnames(lsi)){
	for(c in cov_test_numerical){
		cor_res <- cor.test(x=lsi[[l]], y=cov[[c]])
		p <- cor_res$p.value
		coef <- cor_res$estimate
		batch_df <- rbind(batch_df, list(l,c,coef,p,"continuous"))	
	}
}

colnames(batch_df) <- c("LSIcomponent", "Covariate", "coefficient", "p-value", "class")
batch_df$Covariate <- as.factor(batch_df$Covariate)
batch_df$LSIcomponent <- as.factor(batch_df$LSIcomponent)

# 4. Plot Heatmap of coefficients and p-values

ggplot(batch_df, aes(x=Covariate, y=LSIcomponent, fill=coefficient)) +
	geom_tile() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ArchROutput/Plots/BatchAnalysis_Coefficients.pdf", width = 12, height=10)

ggplot(batch_df, aes(x=Covariate, y=LSIcomponent, fill=-log10(`p-value`))) +
        geom_tile() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("ArchROutput/Plots/BatchAnalysis_p-value.pdf", width = 12, height=10)

# 5. Save project

saveArchRProject(
  ArchRProj = projPM,
  outputDirectory = "ArchROutput/",
)
