library('data.table')

cis_files = c("UKB_460K.blood_PLATELET_COUNT.dat","UKB_460K.blood_PLATELET_COUNT.dat.MHC",
	      "UKB_460K.blood_RBC_DISTRIB_WIDTH.dat","UKB_460K.blood_RBC_DISTRIB_WIDTH.dat.MHC",
	      "UKB_460K.blood_RED_COUNT.dat","UKB_460K.blood_RED_COUNT.dat.MHC",
	      "UKB_460K.blood_WHITE_COUNT.dat","UKB_460K.blood_WHITE_COUNT.dat.MHC",
	      "PASS_Rheumatoid_Arthritis.dat","PASS_Rheumatoid_Arthritis.dat.MHC",
	      "PASS_IBD_deLange2017.dat","PASS_IBD_deLange2017.dat.MHC",
	      "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.dat","UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.dat.MHC",
	      "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.dat","UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.dat.MHC")
trans_files = c("UKB_460K.blood_PLATELET_COUNT.dat","UKB_460K.blood_PLATELET_COUNT.dat.MHC",
              "UKB_460K.blood_RBC_DISTRIB_WIDTH.dat","UKB_460K.blood_RBC_DISTRIB_WIDTH.dat.MHC",
              "UKB_460K.blood_RED_COUNT.dat","UKB_460K.blood_RED_COUNT.dat.MHC",
              "UKB_460K.blood_WHITE_COUNT.dat","UKB_460K.blood_WHITE_COUNT.dat.MHC",
              "PASS_Rheumatoid_Arthritis.dat","PASS_Rheumatoid_Arthritis.dat.MHC",
              "PASS_IBD_deLange2017.dat","PASS_IBD_deLange2017.dat.MHC",
              "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.dat","UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED.dat.MHC",
              "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.dat","UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.dat.MHC")


trans_results =c()
cis_results = c()

cis_dir = "TWAS_association_results/Whole_Blood/cis/"
trans_dir = "TWAS_association_results/Whole_Blood/cis_MatrixeQTL_GBAT_transPCO/"


for (file in trans_files) {
	trans = fread(paste0(trans_dir,file),header = T)
	cis = fread(paste0(cis_dir,file),header = T)
	trans_results = rbind(trans_results,trans)
	cis_results = rbind(cis_results,cis)

	cis$TWAS.P[which(is.na(cis$TWAS.P))] = 1
	trans$TWAS.P[which(is.na(trans$TWAS.P))] = 1
	
        adjusted_cis = p.adjust(cis$TWAS.P,"bonferroni")
	adjusted_trans = p.adjust(trans$TWAS.P,"bonferroni")

	print(file)
	print(sum(adjusted_cis < 0.05))
	print(sum(adjusted_trans < 0.05))
	#print(paste(file,sum(adjusted_trans < 0.05 & adjusted_cis > 0.05)))

	#print(paste(file,sum(adjusted_cis < 0.05 & adjusted_trans > 0.05)))

}
q()
cis_results$TWAS.Z[which(is.na(cis_results$TWAS.Z))] = 0
trans_results$TWAS.Z[which(is.na(trans_results$TWAS.Z))] = 0

cis_results$TWAS.P[which(is.na(cis_results$TWAS.P))] = 1
trans_results$TWAS.P[which(is.na(trans_results$TWAS.P))] = 1

merge_results = merge(cis_results,trans_results, by = "ID")
print(mean(abs(cis_results$TWAS.Z),na.rm =T))
print(mean(abs(trans_results$TWAS.Z),na.rm = T))
print(mean(abs(merge_results$TWAS.Z.x),na.rm = T))
print(mean(abs(merge_results$TWAS.Z.y),na.rm = T))

q()
model_stats = fread("../data/eqtlgen_change_in_r2.txt",header = T)
trans_eqtls = c()
for (row in 1:nrow(eqtlgen_results)) {
	gene = eqtlgen_results$ID[row]
	match = grep(gene,model_stats$gene)
	trans_eqtls = rbind(trans_eqtls,model_stats$'trans eqtls'[match])

}
eqtlgen_results = cbind(eqtlgen_results,trans_eqtls)
print(nrow(cis_results))
print(nrow(eqtlgen_results))

print(mean(abs(eqtlgen_results$TWAS.Z[which(eqtlgen_results$V1 < 1)]),na.rm = T))
print(mean(abs(eqtlgen_results$TWAS.Z[which(eqtlgen_results$V1 < 2)]),na.rm = T))
print(mean(abs(eqtlgen_results$TWAS.Z[which(eqtlgen_results$V1 < 3)]),na.rm = T))
print(mean(abs(eqtlgen_results$TWAS.Z[which(eqtlgen_results$V1 < 4)]),na.rm = T))
print(mean(abs(eqtlgen_results$TWAS.Z[which(eqtlgen_results$V1 > 3)]),na.rm = T))

