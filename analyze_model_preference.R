library(data.table)


results = fread("results_sumstats/Whole_Blood/cis_MatrixeQTL_1e-06_0_mismatches.txt", header = T)
nrow_blup = sum(results$model == 'blup')
nrow_lasso = sum(results$model == 'lasso')
nrow_enet = sum(results$model == 'enet')

num_best_model = c(nrow_blup,nrow_lasso,nrow_enet)
print(num_best_model)
blup_r2 = results$r2[results$model == 'blup']
lasso_r2 = results$r2[results$model == 'lasso']
enet_r2 = results$r2[results$model == 'enet']
rsq_data = list("BLUP" = blup_r2,"lasso" =  lasso_r2,"elastic net" =  enet_r2)

pdf("Count best model.pdf")
barplot(num_best_model, 
        names.arg = c("BLUP", "lasso", "elastic net"), 
        col = "darkgreen", 
        main = "Number of times best model", 
        ylab = "Count", 
        xlab = "models")
dev.off()


pdf("Boxplot of R2 for best models.pdf")
boxplot(rsq_data,
        main = "R-squared Values by best model",
        ylab = "R-squared",
        col = "lightblue")
dev.off()
