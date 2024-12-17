library(data.table)
library('scales')
library(ggplot2)

cis_results = fread("results_sumstats/Whole_Blood/cis_results.txt",header = T)

MatrixeQTL_results = fread("../trans_adapt/results_sumstats/Whole_Blood/cis_MatrixeQTL_1e-06_redo.txt", header = T)

transPCO_results = fread("../transPCO/results_sumstats/Whole_Blood/cis_transPCO_results.txt", header = T)

GBAT_results = fread("results_sumstats/Whole_Blood/cis_GBAT_1e-04.txt",header = T)

cis_results$gene = sapply(strsplit(cis_results$gene, "\\."), function(x) x[1])
GBAT_results$gene = sapply(strsplit(GBAT_results$gene, "\\."), function(x) x[1])

merge_one = merge(cis_results, MatrixeQTL_results,by = 'gene', all.x = T, all.y = T,suffixes = c('.cis','.matrix'))
merge_two = merge(merge_one, transPCO_results,by = 'gene', all.x = T, all.y = T)
merge_three = merge(merge_two,GBAT_results,by = 'gene',all.x = T, all.y = T, suffixes = c('.PCO','.GBAT'))


print(merge_three)

merge_three[is.na(merge_three)] = 0

r2_columns = c("r2.cis", "r2.matrix", "r2.PCO", "r2.GBAT")

# Find the model with the best r2 for each row
merge_three$best_model = apply(merge_three[, ..r2_columns], 1, function(row) {
  # Use which.max to find the index of the maximum r2 value
  colnames(merge_three[, ..r2_columns])[which.max(row)]
})

# Store the best r2 value in another column
merge_three$best_r2 = apply(merge_three[, ..r2_columns], 1, max)

# Print the final result
#print(merge_three[,..r2_columns])
print(sum(merge_three$best_model == 'r2.matrix'))
print(sum(merge_three$best_model == 'r2.PCO'))

best_model_counts = table(merge_three$best_model)
other_models_count <- sum(merge_three$best_model %in% c("r2.matrix", "r2.PCO", "r2.GBAT"))

best_model_counts <- c(best_model_counts, other_models_count)
print(best_model_counts)

model_names = c("Cis", "GBAT", "MatrixeQTL", "transPCO","EGRET")

bar_colors <- c("skyblue", "lightgreen", "orange", "pink", "purple")
bar_colors = c(
  "#377eb8",  # muted blue
  "#e41a1c",  # muted red
  "#4daf4a",  # muted green
  "#984ea3",  # muted purple
  "#ff7f00"  # muted orange
  )
#bar_colors <- c("#FFB3B3", "#B3E6B3", "#B3D1FF", "#FFFFB3", "#CCCCFF")
# Step 2: Create a barplot
pdf("Best Gene Expression Model.pdf")
barplot(best_model_counts,
        main = "Number of Times Each Model is the Best",
        xlab = "Model",
        ylab = "Count",
        col = bar_colors,
        border = "blue",
	names.arg = model_names,
        las = 1)  # las=1 makes the x labels horizontal
dev.off()

merge_three$best_model_name = model_names[merge_three$best_model]  # Label the best model with names
print(merge_three$best_model_names)
# Step 3: Create the scatter plot with ggplot2
pdf("Scatter of r2 by best model.pdf")
ggplot(merge_three, aes(x = r2.cis, y = best_r2, color = best_model)) +
  geom_point() +                                  # Plot points
  scale_color_manual(values = bar_colors,
                     labels = c("Cis", "GBAT", "MatrixeQTL", "transPCO")) +  # Set custom colors for models
  labs(x = "R² Cis", y = "Best R²", color = "Best Model") +   # Labels for axes and legend
  theme_minimal() +                               # Minimal theme for cleaner look
  ggtitle("Scatter plot of R² Cis vs Best Model R²") + 
  theme(axis.title.x = element_text(size = 14),  # Adjust x-axis title size
        axis.title.y = element_text(size = 14),  # Adjust y-axis title size
        axis.text.x = element_text(size = 12),    # Adjust x-axis text size
        axis.text.y = element_text(size = 12),    # Adjust y-axis text size
        legend.title = element_text(size = 14),    # Adjust legend title size
        legend.text = element_text(size = 12),aspect.ratio = 1,text = element_text(size = 14))     # Adjust legend text size

dev.off()
merge_three$r2.matrix[merge_three$r2.matrix < merge_three$r2.cis] = merge_three$r2.cis[merge_three$r2.matrix < merge_three$r2.cis]

merge_three$r2.GBAT[merge_three$r2.GBAT < merge_three$r2.cis] = merge_three$r2.cis[merge_three$r2.GBAT < merge_three$r2.cis]

merge_three$r2.PCO[merge_three$r2.PCO < merge_three$r2.cis] = merge_three$r2.cis[merge_three$r2.PCO < merge_three$r2.cis]

long_r2 <- melt(merge_three[, .(r2.cis, r2.matrix, r2.PCO, r2.GBAT, best_r2)], 
                measure.vars = c("r2.cis", "r2.matrix", "r2.PCO", "r2.GBAT", "best_r2"), 
                variable.name = "Model", value.name = "R2")



pdf("Boxplot_of_R2_Distributions_with_BestR2.pdf")
ggplot(long_r2, aes(x = Model, y = R2, fill = Model)) +
  geom_boxplot() +
  labs(title = "Distribution of R² Values by Model", x = "Model", y = "R²") +
  scale_fill_manual(values = c("r2.cis" = "skyblue", "r2.matrix" = "lightgreen", "r2.PCO" = "orange", "r2.GBAT" = "pink", "best_r2" = "purple")) + 
  scale_x_discrete(labels = c("r2.cis" = "Cis", "r2.matrix" = "MatrixeQTL", "r2.PCO" = "transPCO", "r2.GBAT" = "GBAT", "best_r2" = "Best R²")) +  # Rename x-axis labels
  coord_cartesian(ylim = c(0, 0.2)) +  # Adjust y-axis limits
  theme_minimal() +
  theme(text = element_text(size = 14),aspect.ratio = 1)
dev.off()


sem <- function(x) {
  sd(x) / sqrt(length(x))
}

# Calculate means and SEMs for each model
means <- c(
  mean(merge_three$r2.cis, na.rm = TRUE),
  mean(merge_three$r2.GBAT, na.rm = TRUE),
  mean(merge_three$r2.matrix, na.rm = TRUE),
  mean(merge_three$r2.PCO, na.rm = TRUE),
  mean(merge_three$best_r2, na.rm = TRUE)
)

sems <- c(
  sem(merge_three$r2.cis),
  sem(merge_three$r2.GBAT),
  sem(merge_three$r2.matrix),
  sem(merge_three$r2.PCO),
  sem(merge_three$best_r2)
)

# Model names and colors for consistent labeling
model_names <- c("Cis", "GBAT", "MatrixeQTL", "transPCO", "EGRET")

# Plot means with error bars
pdf("Mean_R2_with_SEM_baseR.pdf")
plot(1:5, means, ylim = c(0.05, 0.065), pch = 16, col = bar_colors,
     xaxt = "n", yaxt = "n",xlab = "Model", ylab = "Mean R²",
     main = "Mean R² for Each Model with SEM")

# Add model names on the x-axis
axis(1, at = 1:5, labels = model_names)
axis(2, at = seq(0.05, 0.065, by = 0.001), labels = seq(0.05, 0.065, by = 0.001))
# Add error bars for SEM
arrows(1:5, means - sems, 1:5, means + sems, angle = 90, code = 3, length = 0.1, col = bar_colors)

dev.off()


print(t.test(merge_three$best_r2, merge_three$r2.cis, paired = FALSE, alternative = "greater")$p.value)
print(t.test(merge_three$best_r2, merge_three$r2.PCO, paired = FALSE, alternative = "greater")$p.value)
print(t.test(merge_three$best_r2, merge_three$r2.GBAT, paired = FALSE, alternative = "greater")$p.value)
print(t.test(merge_three$best_r2, merge_three$r2.matrix, paired = FALSE, alternative = "greater")$p.value)


