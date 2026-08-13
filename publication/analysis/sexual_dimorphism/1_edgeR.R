# This script was developed and run by @blam
# Environment details are broadly missing

library(Seurat)
library(schard)
library(edgeR)

rundate="260629"

sink(paste0("run.",rundate,".log"))
setwd("~/maps/pb")

pb <- h5ad2seurat("../HypoMap3D_full_pb_sub_l4_prelim.h5ad")
pb <- SetIdent(pb, value = pb@meta.data$cell_type)
clusters <- levels(as.factor(pb@meta.data$cell_type))
nsamples <- list()

norm.cutoff <- 0.25
min.lib.size <- 30000

for (j in clusters) {
  print(paste0("working on ", j))
  
  sub <- subset(pb, idents = j)
  
  y <- DGEList(as.matrix(sub@assays$RNA@counts))
  y$samples$meta <- sub@meta.data
  y$samples$group <- as.factor(y$samples$meta$sex_per_sample_inferred)
  
  print("control and M/F only")
  y <- y[, y$samples$meta$IsControl == TRUE & y$samples$meta$sex_per_sample_inferred != "mixed"]
  print(dim(y))
  
  y <- y[, y$samples$lib.size > min.lib.size]
  print("after filtering for min lib size")
  print(dim(y))
  
  y <- calcNormFactors(y)
  y <- y[, abs(y$samples$norm.factors - 1) < norm.cutoff]
  print("after taking out norm outliers")
  print(dim(y))
  
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.size = F]
  y <- calcNormFactors(y)
  
  nsamples[[j]] <- table(y$samples$meta$sex_per_sample_inferred)
  n_datasets <- length(levels(as.factor(y$samples$meta$Dataset)))
  n_sexes <- length(levels(as.factor(y$samples$meta$sex_per_sample_inferred)))
  too_few_replicates <- any(table(y$samples$meta$sex_per_sample_inferred) < 2)
  
  model_case <- if (n_sexes <= 1 | too_few_replicates) {
    "skip"
  } else if (n_datasets > 1) {
    design_matrix <- model.matrix(~y$samples$meta$Dataset + y$samples$meta$sex_per_sample_inferred)
    if (qr(design_matrix)$rank == ncol(design_matrix)) "multivariate" else "univariate_not_full_rank"
  } else {
    "univariate_one_dataset"
  }
  
  switch(model_case,
         skip = {
           print("Skipping pb analysis (only one sex group or not enough replicates)")
         },
         multivariate = {
           print("Full rank, attempting to run multivariate model ~dataset + sex")
           y <- estimateDisp(y, design_matrix, robust = T)
           fit <- glmQLFit(y, design_matrix)
           test <- glmQLFTest(fit, coef = ncol(design_matrix))
           print(topTags(test))
           print(table(decideTests(test)))
           write.table(topTags(test, n = "inf"), file = paste0(j, "_", rundate, "_multivariate_results.txt"), sep = "\t", col.names = NA)
         },
         univariate_not_full_rank = {
           print("The dataset is not full rank, trying univariate model")
           design_matrix <- model.matrix(~y$samples$meta$sex_per_sample_inferred)
           y <- estimateDisp(y, design_matrix, robust = T)
           fit <- glmQLFit(y, design_matrix)
           test <- glmQLFTest(fit, coef = ncol(design_matrix))
           print(topTags(test))
           print(table(decideTests(test)))
           write.table(topTags(test, n = "inf"), file = paste0(j, "_", rundate, "_univariate_non_full_rank_results.txt"), sep = "\t", col.names = NA)
         },
         univariate_one_dataset = {
           print("Cluster only has one dataset, running univariate model")
           design_matrix <- model.matrix(~y$samples$meta$sex_per_sample_inferred)
           y <- estimateDisp(y, design_matrix, robust = T)
           fit <- glmQLFit(y, design_matrix)
           test <- glmQLFTest(fit, coef = ncol(design_matrix))
           print(topTags(test))
           print(table(decideTests(test)))
           write.table(topTags(test, n = "inf"), file = paste0(j, "_", rundate, "_univariate_only_one_dataset_results.txt"), sep = "\t", col.names = NA)
         }
  )
}
sink()
