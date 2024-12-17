library(data.table)

sig_results = data.table(module = numeric(), rsIDs = character(), pval = numeric())
pval_threshold = 0.00001

for (fold in 0:5) {
	num_modules = length(list.files(paste0("transPCO/Whole_Blood/fold_",fold,"/modules/")))
	bed_dir = paste0("transPCO/Whole_Blood/fold_",fold,"/bed_files/")
	dir.create(bed_dir,recursive = T)
	for (module in 1:num_modules) {
		print(module)
		for (chr in 1:22) {
			if (!file.exists(paste0("transPCO/Whole_Blood/fold_",fold,"/PCO_association_results/module_",module,"_chr_",chr,".txt.gz"))) {
				next
			}
			results = fread(paste0("transPCO/Whole_Blood/fold_",fold,"/PCO_association_results/module_",module,"_chr_",chr,".txt.gz"), header = T)
			sig_rows = results[results$pval < pval_threshold,]
			if (nrow(sig_rows) > 0) {
				sig_results = rbind(sig_results, data.table(module = module, sig_rows))
			}
		}

		genes_in_module = fread(paste0("transPCO/Whole_Blood/fold_",fold,"/modules/module_", module, "_gene_list.txt"), header = F)
                genes_in_module = unlist(genes_in_module$V1)

		sig_snps_per_module = sig_results$rsIDs[sig_results$module == module]

                for (gene in genes_in_module) {
                        fwrite(as.matrix(sig_snps_per_module),paste0("transPCO/Whole_Blood/fold_",fold,"/bed_files/",gene,".txt"),col.names = F, row.names = F, sep = '\t',quote = F)
                }

	}

}
