Available_Outcomes<-available_outcomes(access_token = ieugwasr::check_access_token())
ALDO_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/Aldosterone/SummaryStatisticsALDO/ALDO_2SampleMR.txt")
###We imported the above dataset to add new variable called "Phenotype"

#01. ALDO
ALDO_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/Aldosterone/SummaryStatisticsALDO/ALDO_2SampleMR.txt")
ALDO_2SampleMR$Phenotype<-rep("ALDO",length(ALDO_2SampleMR$SNP))
write.table(ALDO_2SampleMR,file = "ALDO_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#02 ANDRO
ANDRO_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/Androstenedione/SummaryStatisticsANDRO/ANDRO_2SampleMR.txt")
ANDRO_2SampleMR$Phenotype<-rep("ANDRO",length(ANDRO_2SampleMR$SNP))
write.table(ANDRO_2SampleMR,file = "ANDRO_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#03 TESTo
TESTO_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/Testosterone/SummaryStatisticsTESTO/TESTO_2SampleMR.txt")
TESTO_2SampleMR$Phenotype<-rep("TESTO",length(TESTO_2SampleMR$SNP))
write.table(TESTO_2SampleMR,file = "TESTO_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#04 ESTRO
ESTRO_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/Estradiol/ESTRO_2SampleMR.txt")
ESTRO_2SampleMR$Phenotype<-rep("ESTRO",length(ESTRO_2SampleMR$SNP))
write.table(ESTRO_2SampleMR,file = "ESTRO_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#05 PROG
PROG_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/Progesterone/SummaryStatisticsPROG/PROG_2SampleMR.txt")
PROG_2SampleMR$Phenotype<-rep("PROG",length(PROG_2SampleMR$SNP))
write.table(PROG_2SampleMR,file = "PROG_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#06 CORT
CORT_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/Cortisol/SummaryStatisticsCORT/CORT_2SampleMR.txt")
CORT_2SampleMR$Phenotype<-rep("CORT",length(CORT_2SampleMR$SNP))
write.table(CORT_2SampleMR,file = "CORT_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#07 OHP17
OHP17_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/17_OH_Progesterone/SummaryStatisticsOHP17/OHP17_2SampleMR.txt")
OHP17_2SampleMR$Phenotype<-rep("OHP17",length(OHP17_2SampleMR$SNP))
write.table(OHP17_2SampleMR,file = "OHP17_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)

#08 DHEAS
DHEAS_2SampleMR <- read.delim("Z:/GENETICS/HyprColoc/DHEAS/SummaryStatisticsDHEAS/DHEAS_2SampleMR.txt")
DHEAS_2SampleMR$Phenotype<-rep("DHEAS",length(DHEAS_2SampleMR$SNP))
write.table(DHEAS_2SampleMR,file = "DHEAS_2SampleMR_updated.txt",sep = "\t",quote = FALSE,row.names = FALSE)
##########################################################################################################################################################################
####creating read_exposure_data files
library(TwoSampleMR)
ALDO<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/ALDO_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")
ANDRO<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/ANDRO_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")
TESTO<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/TESTO_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")
ESTRO<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/ESTRO_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")
PROG<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/PROG_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")
CORT<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/CORT_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")
OHP17<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/OHP17_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")
DHEAS<-read_exposure_data(filename = "Z:/GENETICS/MR_TYPES/DHEAS_2SampleMR_updated.txt",min_pval = 1e-200,clump = FALSE,sep = "\t",phenotype_col = "Phenotype")

###creating extract_outcome_data files
snps_list<-Reduce(union,list(ALDO$SNP,ANDRO$SNP,TESTO$SNP,ESTRO$SNP,PROG$SNP,CORT$SNP,DHEAS$SNP,OHP17$SNP))

PRAD_GWAS<-extract_outcome_data(snps = snps_list,outcomes = as.array("ieu-b-85"))
BRCA_GWAS<-extract_outcome_data(snps = snps_list,outcomes = as.array("ieu-a-1126"))
UCEC_GWAS<-extract_outcome_data(snps = snps_list,outcomes = as.array("ebi-a-GCST006464"))
OV_GWAS<-extract_outcome_data(snps = snps_list,outcomes = as.array("ieu-a-1120"))
COLCA_GWAS_withheader <- read.delim("Z:/GENETICS/SMR-HEIDI/GWAS/CoCa/COLCA_GWAS_withheader.txt")
names(COLCA_GWAS)[names(COLCA_GWAS)=="Phenotype"]<-"outcome"
COLCA_GWAS$id.outcome<-rep("ABCD",length(COLCA_GWAS$SNP))
COLCA_GWAS_2<-subset(COLCA_GWAS,COLCA_GWAS$pval.outcome<0.05)
library(dplyr)
COLCA_GWAS_3<-distinct(.data = COLCA_GWAS,SNP,.keep_all = TRUE)
COLCA_GWAS_4<-read.delim("Z:/GENETICS/MR_TYPES/COLCA_GWAS.txt")

###We had to trnasfer individual exposure and outcome files into QUT-HPC as the HPC job was killed when we load ~21 GB at once.
ALDO_new<-subset(ALDO,ALDO$pval.exposure<5e-6) ##I had to use less stringent p-values for ALDo as the current value did not provide any data
write.table(ALDO_new,file = "ALDO.txt",sep = "\t",quote = FALSE,row.names = FALSE)

ANDRO_new<-subset(ANDRO,ANDRO$pval.exposure<5e-8)
write.table(ANDRO_new,file = "ANDRO.txt",sep = "\t",quote = FALSE,row.names = FALSE)

ESTRO_new<-subset(ESTRO,ESTRO$pval.exposure<5e-8)
write.table(ESTRO_new,file = "ESTRO.txt",sep = "\t",quote = FALSE,row.names = FALSE)

PROG_new<-subset(PROG,PROG$pval.exposure<5e-8)
write.table(PROG_new,file = "PROG.txt",sep = "\t",quote = FALSE,row.names = FALSE)

TESTO_new<-subset(TESTO,TESTO$pval.exposure<5e-8)
write.table(TESTO_new,file = "TESTO.txt",sep = "\t",quote = FALSE,row.names = FALSE)

DHEAS_new<-subset(DHEAS,DHEAS$pval.exposure<5e-8)
write.table(DHEAS_new,file = "DHEAS.txt",sep = "\t",quote = FALSE,row.names = FALSE)

OHP17_new<-subset(OHP17,OHP17$pval.exposure<5e-8)
write.table(OHP17_new,file = "OHP17.txt",sep = "\t",quote = FALSE,row.names = FALSE)

CORT_new<-subset(CORT,CORT$pval.exposure<5e-8)
write.table(CORT_new,file = "CORT.txt",sep = "\t",quote = FALSE,row.names = FALSE)

####Filter GWAS data for p-value<0.05
GWAS_available<-available_outcomes()
View(GWAS_available)
PRAD_new<-subset(Outcome_GWAS_ieu_b_85,Outcome_GWAS_ieu_b_85$pval.outcome<0.05)
BRCA_new<-subset(Outcome_GWAS_ieu_a_1126,Outcome_GWAS_ieu_a_1126$pval.outcome<0.05)
UCEC_new<-subset(Outcome_GWAS_ebi_a,Outcome_GWAS_ebi_a$pval.outcome<0.05)
OV_new<-subset(Outcome_GWAS_ieu_a_1120,Outcome_GWAS_ieu_a_1120$pval.outcome<0.05)
COLCA_GWAS <- read.csv("Z:/GENETICS/SMR-HEIDI/GWAS/CoCa/COLCA_GWAS.txt", sep="")
colnames(COLCA_GWAS)<-c("SNP","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome","outcome","pval.outcome")
COLCA_new<-subset(COLCA_GWAS,COLCA_GWAS$pval.outcome<0.05)
COLCA_GWAS$outcome<-rep("Colorectal",length(COLCA_GWAS$SNP))

write.table(PRAD_new,file = "PRAD_new.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(BRCA_new,file = "BRCA_new.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(UCEC_new,file = "UCEC_new.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(Ov_new,file = "OV_new.txt",sep = "\t",quote=FALSE,row.names=FALSE)
write.table(COLCA_new,file="COLCA_new.txt",sep="\t",quote=FALSE,row.names=FALSE)
####################################################################Prostate cancer#############################################################################
###Change this code for eight hormones:
PRAD_ANDRO<-harmonise_data(ANDRO_new,PRAD_GWAS,action = 2)
PRAD_TESTO<-harmonise_data(TESTO_new,PRAD_GWAS,action = 2)
PRAD_ESTRO<-harmonise_data(ESTRO_new,PRAD_GWAS,action = 2)
PRAD_PROG<-harmonise_data(PROG_new,PRAD_GWAS,action = 2)
PRAD_ALDO<-harmonise_data(ALDO_new,PRAD_GWAS,action = 2)
PRAD_CORT<-harmonise_data(CORT_new,PRAD_GWAS,action = 2)
PRAD_OHP17<-harmonise_data(OHP17_new,PRAD_GWAS,action = 2)
PRAD_DHEAS<-harmonise_data(DHEAS_new,PRAD_GWAS,action = 2)

PRAD_ANDRO_MR_Results<-mr(PRAD_ANDRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_TESTO_MR_Results<-mr(PRAD_TESTO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_ESTRO_MR_Results<-mr(PRAD_ESTRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_PROG_MR_Results<-mr(PRAD_PROG,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_CORT_MR_Results<-mr(PRAD_CORT,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_OHP17_MR_Results<-mr(PRAD_OHP17,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_DHEAS_MR_Results<-mr(PRAD_DHEAS,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_ALDO_MR_Results<-mr(PRAD_ALDO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))


Q_stat_PRAD_ANDRO<-mr_rucker(PRAD_ANDRO)
MR_PRESSO_PRAD_ANDRO<-run_mr_presso(dat = PRAD_ANDRO)
MR_egger_PRAD_ANDRO<-mr_egger_regression(b_exp = PRAD_ANDRO$beta.exposure,b_out = PRAD_ANDRO$beta.outcome,se_exp = PRAD_ANDRO$se.exposure,se_out = PRAD_ANDRO$se.outcome)

Q_stat_PRAD_TESTO<-mr_rucker(PRAD_TESTO)
MR_PRESSO_PRAD_TESTO<-run_mr_presso(dat = PRAD_TESTO)
MR_egger_PRAD_TESTO<-mr_egger_regression(b_exp = PRAD_TESTO$beta.exposure,b_out = PRAD_TESTO$beta.outcome,se_exp = PRAD_TESTO$se.exposure,se_out = PRAD_TESTO$se.outcome)

Q_stat_PRAD_ESTRO<-mr_rucker(PRAD_ESTRO)
MR_PRESSO_PRAD_ESTRO<-run_mr_presso(dat = PRAD_ESTRO)
MR_egger_PRAD_ESTRO<-mr_egger_regression(b_exp = PRAD_ESTRO$beta.exposure,b_out = PRAD_ESTRO$beta.outcome,se_exp = PRAD_ESTRO$se.exposure,se_out = PRAD_ESTRO$se.outcome)

Q_stat_PRAD_PROG<-mr_rucker(PRAD_PROG)
MR_PRESSO_PRAD_PROG<-run_mr_presso(dat = PRAD_PROG)
MR_egger_PRAD_PROG<-mr_egger_regression(b_exp = PRAD_PROG$beta.exposure,b_out = PRAD_PROG$beta.outcome,se_exp = PRAD_PROG$se.exposure,se_out = PRAD_PROG$se.outcome)

Q_stat_PRAD_CORT<-mr_rucker(PRAD_CORT)
MR_PRESSO_PRAD_CORT<-run_mr_presso(dat = PRAD_CORT)
MR_egger_PRAD_CORT<-mr_egger_regression(b_exp = PRAD_CORT$beta.exposure,b_out = PRAD_CORT$beta.outcome,se_exp = PRAD_CORT$se.exposure,se_out = PRAD_CORT$se.outcome)

Q_stat_PRAD_OHP17<-mr_rucker(PRAD_OHP17)
MR_PRESSO_PRAD_OHP17<-run_mr_presso(dat = PRAD_OHP17)
MR_egger_PRAD_OHP17<-mr_egger_regression(b_exp = PRAD_OHP17$beta.exposure,b_out = PRAD_OHP17$beta.outcome,se_exp = PRAD_OHP17$se.exposure,se_out = PRAD_OHP17$se.outcome)

Q_stat_PRAD_DHEAS<-mr_rucker(PRAD_DHEAS)
MR_PRESSO_PRAD_DHEAS<-run_mr_presso(dat = PRAD_DHEAS)
MR_egger_PRAD_DHEAS<-mr_egger_regression(b_exp = PRAD_DHEAS$beta.exposure,b_out = PRAD_DHEAS$beta.outcome,se_exp = PRAD_DHEAS$se.exposure,se_out = PRAD_DHEAS$se.outcome)

Q_stat_PRAD_ALDO<-mr_rucker(PRAD_ALDO)
MR_PRESSO_PRAD_ALDO<-run_mr_presso(dat = PRAD_ALDO)
MR_egger_PRAD_ALDO<-mr_egger_regression(b_exp = PRAD_ALDO$beta.exposure,b_out = PRAD_ALDO$beta.outcome,se_exp = PRAD_ALDO$se.exposure,se_out = PRAD_ALDO$se.outcome)


write.table(PRAD_ANDRO_MR_Results,file="PRAD_ANDRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(PRAD_ALDO_MR_Results,file="PRAD_ALDO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(PRAD_ESTRO_MR_Results,file="PRAD_ESTRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(PRAD_PROG_MR_Results,file="PRAD_PROG_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(PRAD_TESTO_MR_Results,file="PRAD_TESTO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(PRAD_CORT_MR_Results,file="PRAD_CORT_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(PRAD_OHP17_MR_Results,file="PRAD_OHP17_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(PRAD_DHEAS_MR_Results,file="PRAD_DHEAS_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)

#for bidirectional analysis

###Formatting data
PRAD_GWAS_2<-PRAD_GWAS
names(PRAD_GWAS_2)<-gsub(x = names(PRAD_GWAS_2), pattern = "outcome", replacement = "exposure")
BRCA_GWAS_2<-BRCA_GWAS
names(BRCA_GWAS_2)<-gsub(x = names(BRCA_GWAS_2), pattern = "outcome", replacement = "exposure")
UCEC_GWAS_2<-UCEC_GWAS
names(UCEC_GWAS_2)<-gsub(x = names(UCEC_GWAS_2), pattern = "outcome", replacement = "exposure")
OV_GWAS_2<-OV_GWAS
names(OV_GWAS_2)<-gsub(x = names(OV_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")

ESTRO_new_2<-ESTRO_new
names(ESTRO_new_2)<-gsub(x = names(ESTRO_new_2), pattern = "exposure", replacement = "outcome")
PROG_new_2<-PROG_new
names(PROG_new_2)<-gsub(x = names(PROG_new_2), pattern = "exposure", replacement = "outcome")
TESTO_new_2<-TESTO_new
names(TESTO_new_2)<-gsub(x = names(TESTO_new_2), pattern = "exposure", replacement = "outcome")
ANDRO_new_2<-ANDRO_new
names(ANDRO_new_2)<-gsub(x = names(ANDRO_new_2), pattern = "exposure", replacement = "outcome")
ALDO_new_2<-ALDO_new
names(ALDO_new_2)<-gsub(x = names(ALDO_new_2), pattern = "exposure", replacement = "outcome")
OHP17_new_2<-OHP17_new
names(OHP17_new_2)<-gsub(x = names(OHP17_new_2), pattern = "exposure", replacement = "outcome")
DHEAS_new_2<-DHEAS_new
names(DHEAS_new_2)<-gsub(x = names(DHEAS_new_2), pattern = "exposure", replacement = "outcome")
CORT_new_2<-CORT_new
names(CORT_new_2)<-gsub(x = names(CORT_new_2), pattern = "exposure", replacement = "outcome")

PRAD_ANDRO_2<-harmonise_data(PRAD_GWAS_2,ANDRO_new_2,action = 2)
PRAD_ANDRO_2<-steiger_filtering(PRAD_ANDRO_2)
PRAD_TESTO_2<-harmonise_data(TESTO_new,PRAD_GWAS,action = 2)
PRAD_ESTRO_2<-harmonise_data(ESTRO_new,PRAD_GWAS,action = 2)
PRAD_PROG_2<-harmonise_data(PROG_new,PRAD_GWAS,action = 2)
PRAD_ALDO_2<-harmonise_data(ALDO_new,PRAD_GWAS,action = 2)
PRAD_ALDO_2<-steiger_filtering(PRAD_ALDO_2)
PRAD_CORT_2<-harmonise_data(CORT_new,PRAD_GWAS,action = 2)
PRAD_OHP17_2<-harmonise_data(OHP17_new,PRAD_GWAS,action = 2)
PRAD_DHEAS_2<-harmonise_data(DHEAS_new,PRAD_GWAS,action = 2)

PRAD_ANDRO_MR_Results2<-mr(PRAD_ANDRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_TESTO_MR_Results2<-mr(PRAD_TESTO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_ESTRO_MR_Results2<-mr(PRAD_ESTRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_PROG_MR_Results2<-mr(PRAD_PROG_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_CORT_MR_Results2<-mr(PRAD_CORT_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_OHP17_MR_Results2<-mr(PRAD_OHP17_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_DHEAS_MR_Results2<-mr(PRAD_DHEAS_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
PRAD_ALDO_MR_Results2<-mr(PRAD_ALDO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))


Q_stat_PRAD_ANDRO2<-mr_rucker(PRAD_ANDRO_2)
MR_PRESSO_PRAD_ANDRO2<-run_mr_presso(dat = PRAD_ANDRO_2)
MR_egger_PRAD_ANDRO2<-mr_egger_regression(b_exp = PRAD_ANDRO_2$beta.exposure,b_out = PRAD_ANDRO_2$beta.outcome,se_exp = PRAD_ANDRO_2$se.exposure,se_out = PRAD_ANDRO_2$se.outcome)

Q_stat_PRAD_TESTO2<-mr_rucker(PRAD_TESTO_2)
MR_PRESSO_PRAD_TESTO2<-run_mr_presso(dat = PRAD_TESTO_2)
MR_egger_PRAD_TESTO2<-mr_egger_regression(b_exp = PRAD_TESTO_2$beta.exposure,b_out = PRAD_TESTO_2$beta.outcome,se_exp = PRAD_TESTO_2$se.exposure,se_out = PRAD_TESTO_2$se.outcome)

Q_stat_PRAD_ESTRO2<-mr_rucker(PRAD_ESTRO_2)
MR_PRESSO_PRAD_ESTRO2<-run_mr_presso(dat = PRAD_ESTRO_2)
MR_egger_PRAD_ESTRO2<-mr_egger_regression(b_exp = PRAD_ESTRO_2$beta.exposure,b_out = PRAD_ESTRO_2$beta.outcome,se_exp = PRAD_ESTRO_2$se.exposure,se_out = PRAD_ESTRO_2$se.outcome)

Q_stat_PRAD_PROG2<-mr_rucker(PRAD_PROG_2)
MR_PRESSO_PRAD_PROG2<-run_mr_presso(dat = PRAD_PROG_2)
MR_egger_PRAD_PROG2<-mr_egger_regression(b_exp = PRAD_PROG_2$beta.exposure,b_out = PRAD_PROG_2$beta.outcome,se_exp = PRAD_PROG_2$se.exposure,se_out = PRAD_PROG_2$se.outcome)

Q_stat_PRAD_CORT2<-mr_rucker(PRAD_CORT_2)
MR_PRESSO_PRAD_CORT2<-run_mr_presso(dat = PRAD_CORT_2)
MR_egger_PRAD_CORT2<-mr_egger_regression(b_exp = PRAD_CORT_2$beta.exposure,b_out = PRAD_CORT_2$beta.outcome,se_exp = PRAD_CORT_2$se.exposure,se_out = PRAD_CORT_2$se.outcome)

Q_stat_PRAD_OHP172<-mr_rucker(PRAD_OHP17_2)
MR_PRESSO_PRAD_OHP172<-run_mr_presso(dat = PRAD_OHP17_2)
MR_egger_PRAD_OHP17<-mr_egger_regression(b_exp = PRAD_OHP17_2$beta.exposure,b_out = PRAD_OHP17_2$beta.outcome,se_exp = PRAD_OHP17_2$se.exposure,se_out = PRAD_OHP17_2$se.outcome)

Q_stat_PRAD_DHEAS2<-mr_rucker(PRAD_DHEAS_2)
MR_PRESSO_PRAD_DHEAS2<-run_mr_presso(dat = PRAD_DHEAS_2)
MR_egger_PRAD_DHEAS2<-mr_egger_regression(b_exp = PRAD_DHEAS_2$beta.exposure,b_out = PRAD_DHEAS_2$beta.outcome,se_exp = PRAD_DHEAS_2$se.exposure,se_out = PRAD_DHEAS_2$se.outcome)

Q_stat_PRAD_ALDO2<-mr_rucker(PRAD_ALDO_2)
MR_PRESSO_PRAD_ALDO2<-run_mr_presso(dat = PRAD_ALDO_2)
MR_egger_PRAD_ALDO2<-mr_egger_regression(b_exp = PRAD_ALDO_2$beta.exposure,b_out = PRAD_ALDO_2$beta.outcome,se_exp = PRAD_ALDO_2$se.exposure,se_out = PRAD_ALDO_2$se.outcome)


##########################################################Breast cancer######################################################################
###Change this code for eight hormones:
BRCA_ANDRO<-harmonise_data(ANDRO_new,BRCA_GWAS,action = 2)
BRCA_TESTO<-harmonise_data(TESTO_new,BRCA_GWAS,action = 2)
BRCA_ESTRO<-harmonise_data(ESTRO_new,BRCA_GWAS,action = 2)
BRCA_PROG<-harmonise_data(PROG_new,BRCA_GWAS,action = 2)
BRCA_ALDO<-harmonise_data(ALDO_new,BRCA_GWAS,action = 2)
BRCA_CORT<-harmonise_data(CORT_new,BRCA_GWAS,action = 2)
BRCA_OHP17<-harmonise_data(OHP17_new,BRCA_GWAS,action = 2)
BRCA_DHEAS<-harmonise_data(DHEAS_new,BRCA_GWAS,action = 2)

####Read above files (using read.delim, see toHPC.pbs file) into HPC and ran following codes
BRCA_ANDRO_MR_Results<-mr(BRCA_ANDRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_TESTO_MR_Results<-mr(BRCA_TESTO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_ESTRO_MR_Results<-mr(BRCA_ESTRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_PROG_MR_Results<-mr(BRCA_PROG,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_CORT_MR_Results<-mr(BRCA_CORT,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_OHP17_MR_Results<-mr(BRCA_OHP17,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_DHEAS_MR_Results<-mr(BRCA_DHEAS,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_ALDO_MR_Results<-mr(BRCA_ALDO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))

Q_stat_BRCA_ANDRO<-mr_rucker(BRCA_ANDRO)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_ANDRO)
MR_egger_BRCA_ANDRO<-mr_egger_regression(b_exp = BRCA_ANDRO$beta.exposure,b_out = BRCA_ANDRO$beta.outcome,se_exp = BRCA_ANDRO$se.exposure,se_out = BRCA_ANDRO$se.outcome)

Q_stat_BRCA_TESTO<-mr_rucker(BRCA_TESTO)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_TESTO)
MR_egger_BRCA_TESTO<-mr_egger_regression(b_exp = BRCA_TESTO$beta.exposure,b_out = BRCA_TESTO$beta.outcome,se_exp = BRCA_TESTO$se.exposure,se_out = BRCA_TESTO$se.outcome)

Q_stat_BRCA_ESTRO<-mr_rucker(BRCA_ESTRO)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_ESTRO)
MR_egger_BRCA_ESTRO<-mr_egger_regression(b_exp = BRCA_ESTRO$beta.exposure,b_out = BRCA_ESTRO$beta.outcome,se_exp = BRCA_ESTRO$se.exposure,se_out = BRCA_ESTRO$se.outcome)

Q_stat_BRCA_PROG<-mr_rucker(BRCA_PROG)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_PROG)
MR_egger_BRCA_PROG<-mr_egger_regression(b_exp = BRCA_PROG$beta.exposure,b_out = BRCA_PROG$beta.outcome,se_exp = BRCA_PROG$se.exposure,se_out = BRCA_PROG$se.outcome)

Q_stat_BRCA_CORT<-mr_rucker(BRCA_CORT)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_CORT)
MR_egger_BRCA_CORT<-mr_egger_regression(b_exp = BRCA_CORT$beta.exposure,b_out = BRCA_CORT$beta.outcome,se_exp = BRCA_CORT$se.exposure,se_out = BRCA_CORT$se.outcome)

Q_stat_BRCA_OHP17<-mr_rucker(BRCA_OHP17)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_OHP17)
MR_egger_BRCA_OHP17<-mr_egger_regression(b_exp = BRCA_OHP17$beta.exposure,b_out = BRCA_OHP17$beta.outcome,se_exp = BRCA_OHP17$se.exposure,se_out = BRCA_OHP17$se.outcome)

Q_stat_BRCA_DHEAS<-mr_rucker(BRCA_DHEAS)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_DHEAS)
MR_egger_BRCA_DHEAS<-mr_egger_regression(b_exp = BRCA_DHEAS$beta.exposure,b_out = BRCA_DHEAS$beta.outcome,se_exp = BRCA_DHEAS$se.exposure,se_out = BRCA_DHEAS$se.outcome)

Q_stat_BRCA_ALDO<-mr_rucker(BRCA_ALDO)
MR_PRESSO_BRCA_AQNDRO<-run_mr_presso(dat = BRCA_ALDO)
MR_egger_BRCA_ALDO<-mr_egger_regression(b_exp = BRCA_ALDO$beta.exposure,b_out = BRCA_ALDO$beta.outcome,se_exp = BRCA_ALDO$se.exposure,se_out = BRCA_ALDO$se.outcome)



write.table(BRCA_ANDRO_MR_Results,file="BRCA_ANDRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(BRCA_ALDO_MR_Results,file="BRCA_ALDO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(BRCA_ESTRO_MR_Results,file="BRCA_ESTRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(BRCA_PROG_MR_Results,file="BRCA_PROG_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(BRCA_TESTO_MR_Results,file="BRCA_TESTO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(BRCA_CORT_MR_Results,file="BRCA_CORT_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(BRCA_OHP17_MR_Results,file="BRCA_OHP17_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(BRCA_DHEAS_MR_Results,file="BRCA_DHEAS_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)

#for bidirectional analysis
BRCA_GWAS_2<-BRCA_GWAS
names(BRCA_GWAS_2)<-gsub(x = names(BRCA_GWAS_2), pattern = "outcome", replacement = "exposure")
BRCA_GWAS_2<-BRCA_GWAS
names(BRCA_GWAS_2)<-gsub(x = names(BRCA_GWAS_2), pattern = "outcome", replacement = "exposure")
UCEC_GWAS_2<-UCEC_GWAS
names(UCEC_GWAS_2)<-gsub(x = names(UCEC_GWAS_2), pattern = "outcome", replacement = "exposure")
OV_GWAS_2<-OV_GWAS
names(OV_GWAS_2)<-gsub(x = names(OV_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")

ESTRO_new_2<-ESTRO_new
names(ESTRO_new_2)<-gsub(x = names(ESTRO_new_2), pattern = "exposure", replacement = "outcome")
PROG_new_2<-PROG_new
names(PROG_new_2)<-gsub(x = names(PROG_new_2), pattern = "exposure", replacement = "outcome")
TESTO_new_2<-TESTO_new
names(TESTO_new_2)<-gsub(x = names(TESTO_new_2), pattern = "exposure", replacement = "outcome")
ANDRO_new_2<-ANDRO_new
names(ANDRO_new_2)<-gsub(x = names(ANDRO_new_2), pattern = "exposure", replacement = "outcome")
ALDO_new_2<-ALDO_new
names(ALDO_new_2)<-gsub(x = names(ALDO_new_2), pattern = "exposure", replacement = "outcome")
OHP17_new_2<-OHP17_new
names(OHP17_new_2)<-gsub(x = names(OHP17_new_2), pattern = "exposure", replacement = "outcome")
DHEAS_new_2<-DHEAS_new
names(DHEAS_new_2)<-gsub(x = names(DHEAS_new_2), pattern = "exposure", replacement = "outcome")
CORT_new_2<-CORT_new
names(CORT_new_2)<-gsub(x = names(CORT_new_2), pattern = "exposure", replacement = "outcome")

BRCA_ANDRO_2<-harmonise_data(BRCA_GWAS_2,ANDRO_new_2,action = 2)
BRCA_TESTO_2<-harmonise_data(TESTO_new,BRCA_GWAS,action = 2)
BRCA_ESTRO_2<-harmonise_data(ESTRO_new,BRCA_GWAS,action = 2)
BRCA_PROG_2<-harmonise_data(PROG_new,BRCA_GWAS,action = 2)
BRCA_ALDO_2<-harmonise_data(ALDO_new,BRCA_GWAS,action = 2)
BRCA_CORT_2<-harmonise_data(CORT_new,BRCA_GWAS,action = 2)
BRCA_OHP17_2<-harmonise_data(OHP17_new,BRCA_GWAS,action = 2)
BRCA_DHEAS_2<-harmonise_data(DHEAS_new,BRCA_GWAS,action = 2)

BRCA_ANDRO_MR_Results2<-mr(BRCA_ANDRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_TESTO_MR_Results2<-mr(BRCA_TESTO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_ESTRO_MR_Results2<-mr(BRCA_ESTRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_PROG_MR_Results2<-mr(BRCA_PROG_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_CORT_MR_Results2<-mr(BRCA_CORT_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_OHP17_MR_Results2<-mr(BRCA_OHP17_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_DHEAS_MR_Results2<-mr(BRCA_DHEAS_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
BRCA_ALDO_MR_Results2<-mr(BRCA_ALDO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))


Q_stat_BRCA_ANDRO2<-mr_rucker(BRCA_ANDRO_2)
MR_PRESSO_BRCA_ANDRO2<-run_mr_presso(dat = BRCA_ANDRO_2)
MR_egger_BRCA_ANDRO2<-mr_egger_regression(b_exp = BRCA_ANDRO_2$beta.exposure,b_out = BRCA_ANDRO_2$beta.outcome,se_exp = BRCA_ANDRO_2$se.exposure,se_out = BRCA_ANDRO_2$se.outcome)

Q_stat_BRCA_TESTO2<-mr_rucker(BRCA_TESTO_2)
MR_PRESSO_BRCA_TESTO2<-run_mr_presso(dat = BRCA_TESTO_2)
MR_egger_BRCA_TESTO2<-mr_egger_regression(b_exp = BRCA_TESTO_2$beta.exposure,b_out = BRCA_TESTO_2$beta.outcome,se_exp = BRCA_TESTO_2$se.exposure,se_out = BRCA_TESTO_2$se.outcome)

Q_stat_BRCA_ESTRO2<-mr_rucker(BRCA_ESTRO_2)
MR_PRESSO_BRCA_ESTRO2<-run_mr_presso(dat = BRCA_ESTRO_2)
MR_egger_BRCA_ESTRO2<-mr_egger_regression(b_exp = BRCA_ESTRO_2$beta.exposure,b_out = BRCA_ESTRO_2$beta.outcome,se_exp = BRCA_ESTRO_2$se.exposure,se_out = BRCA_ESTRO_2$se.outcome)

Q_stat_BRCA_PROG2<-mr_rucker(BRCA_PROG_2)
MR_PRESSO_BRCA_PROG2<-run_mr_presso(dat = BRCA_PROG_2)
MR_egger_BRCA_PROG2<-mr_egger_regression(b_exp = BRCA_PROG_2$beta.exposure,b_out = BRCA_PROG_2$beta.outcome,se_exp = BRCA_PROG_2$se.exposure,se_out = BRCA_PROG_2$se.outcome)

Q_stat_BRCA_CORT2<-mr_rucker(BRCA_CORT_2)
MR_PRESSO_BRCA_CORT2<-run_mr_presso(dat = BRCA_CORT_2)
MR_egger_BRCA_CORT2<-mr_egger_regression(b_exp = BRCA_CORT_2$beta.exposure,b_out = BRCA_CORT_2$beta.outcome,se_exp = BRCA_CORT_2$se.exposure,se_out = BRCA_CORT_2$se.outcome)

Q_stat_BRCA_OHP172<-mr_rucker(BRCA_OHP17_2)
MR_PRESSO_BRCA_OHP172<-run_mr_presso(dat = BRCA_OHP17_2)
MR_egger_BRCA_OHP17<-mr_egger_regression(b_exp = BRCA_OHP17_2$beta.exposure,b_out = BRCA_OHP17_2$beta.outcome,se_exp = BRCA_OHP17_2$se.exposure,se_out = BRCA_OHP17_2$se.outcome)

Q_stat_BRCA_DHEAS2<-mr_rucker(BRCA_DHEAS_2)
MR_PRESSO_BRCA_DHEAS2<-run_mr_presso(dat = BRCA_DHEAS_2)
MR_egger_BRCA_DHEAS2<-mr_egger_regression(b_exp = BRCA_DHEAS_2$beta.exposure,b_out = BRCA_DHEAS_2$beta.outcome,se_exp = BRCA_DHEAS_2$se.exposure,se_out = BRCA_DHEAS_2$se.outcome)

Q_stat_BRCA_ALDO2<-mr_rucker(BRCA_ALDO_2)
MR_PRESSO_BRCA_ALDO2<-run_mr_presso(dat = BRCA_ALDO_2)
MR_egger_BRCA_ALDO2<-mr_egger_regression(b_exp = BRCA_ALDO_2$beta.exposure,b_out = BRCA_ALDO_2$beta.outcome,se_exp = BRCA_ALDO_2$se.exposure,se_out = BRCA_ALDO_2$se.outcome)


###############################################################Ovarian cancer################################################################################
###Change this code for eight hormones:
OV_ANDRO<-harmonise_data(ANDRO_new,OV_GWAS,action = 2)
OV_TESTO<-harmonise_data(TESTO_new,OV_GWAS,action = 2)
OV_ESTRO<-harmonise_data(ESTRO_new,OV_GWAS,action = 2)
OV_PROG<-harmonise_data(PROG_new,OV_GWAS,action = 2)
OV_ALDO<-harmonise_data(ALDO_new,OV_GWAS,action = 2)
OV_CORT<-harmonise_data(CORT_new,OV_GWAS,action = 2)
OV_OHP17<-harmonise_data(OHP17_new,OV_GWAS,action = 2)
OV_DHEAS<-harmonise_data(DHEAS_new,OV_GWAS,action = 2)


OV_ANDRO_MR_Results<-mr(OV_ANDRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_TESTO_MR_Results<-mr(OV_TESTO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_ESTRO_MR_Results<-mr(OV_ESTRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_PROG_MR_Results<-mr(OV_PROG,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_CORT_MR_Results<-mr(OV_CORT,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_OHP17_MR_Results<-mr(OV_OHP17,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_DHEAS_MR_Results<-mr(OV_DHEAS,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_ALDO_MR_Results<-mr(OV_ALDO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))

Q_stat_OV_ANDRO<-mr_rucker(OV_ANDRO)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_ANDRO)
MR_egger_OV_ANDRO<-mr_egger_regression(b_exp = OV_ANDRO$beta.exposure,b_out = OV_ANDRO$beta.outcome,se_exp = OV_ANDRO$se.exposure,se_out = OV_ANDRO$se.outcome)

Q_stat_OV_TESTO<-mr_rucker(OV_TESTO)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_TESTO)
MR_egger_OV_TESTO<-mr_egger_regression(b_exp = OV_TESTO$beta.exposure,b_out = OV_TESTO$beta.outcome,se_exp = OV_TESTO$se.exposure,se_out = OV_TESTO$se.outcome)

Q_stat_OV_ESTRO<-mr_rucker(OV_ESTRO)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_ESTRO)
MR_egger_OV_ESTRO<-mr_egger_regression(b_exp = OV_ESTRO$beta.exposure,b_out = OV_ESTRO$beta.outcome,se_exp = OV_ESTRO$se.exposure,se_out = OV_ESTRO$se.outcome)

Q_stat_OV_PROG<-mr_rucker(OV_PROG)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_PROG)
MR_egger_OV_PROG<-mr_egger_regression(b_exp = OV_PROG$beta.exposure,b_out = OV_PROG$beta.outcome,se_exp = OV_PROG$se.exposure,se_out = OV_PROG$se.outcome)

Q_stat_OV_CORT<-mr_rucker(OV_CORT)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_CORT)
MR_egger_OV_CORT<-mr_egger_regression(b_exp = OV_CORT$beta.exposure,b_out = OV_CORT$beta.outcome,se_exp = OV_CORT$se.exposure,se_out = OV_CORT$se.outcome)

Q_stat_OV_OHP17<-mr_rucker(OV_OHP17)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_OHP17)
MR_egger_OV_OHP17<-mr_egger_regression(b_exp = OV_OHP17$beta.exposure,b_out = OV_OHP17$beta.outcome,se_exp = OV_OHP17$se.exposure,se_out = OV_OHP17$se.outcome)

Q_stat_OV_DHEAS<-mr_rucker(OV_DHEAS)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_DHEAS)
MR_egger_OV_DHEAS<-mr_egger_regression(b_exp = OV_DHEAS$beta.exposure,b_out = OV_DHEAS$beta.outcome,se_exp = OV_DHEAS$se.exposure,se_out = OV_DHEAS$se.outcome)

Q_stat_OV_ALDO<-mr_rucker(OV_ALDO)
MR_PRESSO_OV_AQNDRO<-run_mr_presso(dat = OV_ALDO)
MR_egger_OV_ALDO<-mr_egger_regression(b_exp = OV_ALDO$beta.exposure,b_out = OV_ALDO$beta.outcome,se_exp = OV_ALDO$se.exposure,se_out = OV_ALDO$se.outcome)


write.table(OV_ANDRO_MR_Results,file="OV_ANDRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(OV_ALDO_MR_Results,file="OV_ALDO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(OV_ESTRO_MR_Results,file="OV_ESTRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(OV_PROG_MR_Results,file="OV_PROG_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(OV_TESTO_MR_Results,file="OV_TESTO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(OV_CORT_MR_Results,file="OV_CORT_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(OV_OHP17_MR_Results,file="OV_OHP17_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(OV_DHEAS_MR_Results,file="OV_DHEAS_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)

#for bidirectional analysis
OV_GWAS_2<-OV_GWAS
names(OV_GWAS_2)<-gsub(x = names(OV_GWAS_2), pattern = "outcome", replacement = "exposure")
OV_GWAS_2<-OV_GWAS
names(OV_GWAS_2)<-gsub(x = names(OV_GWAS_2), pattern = "outcome", replacement = "exposure")
UCEC_GWAS_2<-UCEC_GWAS
names(UCEC_GWAS_2)<-gsub(x = names(UCEC_GWAS_2), pattern = "outcome", replacement = "exposure")
OV_GWAS_2<-OV_GWAS
names(OV_GWAS_2)<-gsub(x = names(OV_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")

ESTRO_new_2<-ESTRO_new
names(ESTRO_new_2)<-gsub(x = names(ESTRO_new_2), pattern = "exposure", replacement = "outcome")
PROG_new_2<-PROG_new
names(PROG_new_2)<-gsub(x = names(PROG_new_2), pattern = "exposure", replacement = "outcome")
TESTO_new_2<-TESTO_new
names(TESTO_new_2)<-gsub(x = names(TESTO_new_2), pattern = "exposure", replacement = "outcome")
ANDRO_new_2<-ANDRO_new
names(ANDRO_new_2)<-gsub(x = names(ANDRO_new_2), pattern = "exposure", replacement = "outcome")
ALDO_new_2<-ALDO_new
names(ALDO_new_2)<-gsub(x = names(ALDO_new_2), pattern = "exposure", replacement = "outcome")
OHP17_new_2<-OHP17_new
names(OHP17_new_2)<-gsub(x = names(OHP17_new_2), pattern = "exposure", replacement = "outcome")
DHEAS_new_2<-DHEAS_new
names(DHEAS_new_2)<-gsub(x = names(DHEAS_new_2), pattern = "exposure", replacement = "outcome")
CORT_new_2<-CORT_new
names(CORT_new_2)<-gsub(x = names(CORT_new_2), pattern = "exposure", replacement = "outcome")

OV_ANDRO_2<-harmonise_data(OV_GWAS_2,ANDRO_new_2,action = 2)
OV_TESTO_2<-harmonise_data(TESTO_new,OV_GWAS,action = 2)
OV_ESTRO_2<-harmonise_data(ESTRO_new,OV_GWAS,action = 2)
OV_PROG_2<-harmonise_data(PROG_new,OV_GWAS,action = 2)
OV_ALDO_2<-harmonise_data(ALDO_new,OV_GWAS,action = 2)
OV_CORT_2<-harmonise_data(CORT_new,OV_GWAS,action = 2)
OV_OHP17_2<-harmonise_data(OHP17_new,OV_GWAS,action = 2)
OV_DHEAS_2<-harmonise_data(DHEAS_new,OV_GWAS,action = 2)

OV_ANDRO_MR_Results2<-mr(OV_ANDRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_TESTO_MR_Results2<-mr(OV_TESTO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_ESTRO_MR_Results2<-mr(OV_ESTRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_PROG_MR_Results2<-mr(OV_PROG_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_CORT_MR_Results2<-mr(OV_CORT_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_OHP17_MR_Results2<-mr(OV_OHP17_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_DHEAS_MR_Results2<-mr(OV_DHEAS_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
OV_ALDO_MR_Results2<-mr(OV_ALDO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))


Q_stat_OV_ANDRO2<-mr_rucker(OV_ANDRO_2)
MR_PRESSO_OV_ANDRO2<-run_mr_presso(dat = OV_ANDRO_2)
MR_egger_OV_ANDRO2<-mr_egger_regression(b_exp = OV_ANDRO_2$beta.exposure,b_out = OV_ANDRO_2$beta.outcome,se_exp = OV_ANDRO_2$se.exposure,se_out = OV_ANDRO_2$se.outcome)

Q_stat_OV_TESTO2<-mr_rucker(OV_TESTO_2)
MR_PRESSO_OV_TESTO2<-run_mr_presso(dat = OV_TESTO_2)
MR_egger_OV_TESTO2<-mr_egger_regression(b_exp = OV_TESTO_2$beta.exposure,b_out = OV_TESTO_2$beta.outcome,se_exp = OV_TESTO_2$se.exposure,se_out = OV_TESTO_2$se.outcome)

Q_stat_OV_ESTRO2<-mr_rucker(OV_ESTRO_2)
MR_PRESSO_OV_ESTRO2<-run_mr_presso(dat = OV_ESTRO_2)
MR_egger_OV_ESTRO2<-mr_egger_regression(b_exp = OV_ESTRO_2$beta.exposure,b_out = OV_ESTRO_2$beta.outcome,se_exp = OV_ESTRO_2$se.exposure,se_out = OV_ESTRO_2$se.outcome)

Q_stat_OV_PROG2<-mr_rucker(OV_PROG_2)
MR_PRESSO_OV_PROG2<-run_mr_presso(dat = OV_PROG_2)
MR_egger_OV_PROG2<-mr_egger_regression(b_exp = OV_PROG_2$beta.exposure,b_out = OV_PROG_2$beta.outcome,se_exp = OV_PROG_2$se.exposure,se_out = OV_PROG_2$se.outcome)

Q_stat_OV_CORT2<-mr_rucker(OV_CORT_2)
MR_PRESSO_OV_CORT2<-run_mr_presso(dat = OV_CORT_2)
MR_egger_OV_CORT2<-mr_egger_regression(b_exp = OV_CORT_2$beta.exposure,b_out = OV_CORT_2$beta.outcome,se_exp = OV_CORT_2$se.exposure,se_out = OV_CORT_2$se.outcome)

Q_stat_OV_OHP172<-mr_rucker(OV_OHP17_2)
MR_PRESSO_OV_OHP172<-run_mr_presso(dat = OV_OHP17_2)
MR_egger_OV_OHP17<-mr_egger_regression(b_exp = OV_OHP17_2$beta.exposure,b_out = OV_OHP17_2$beta.outcome,se_exp = OV_OHP17_2$se.exposure,se_out = OV_OHP17_2$se.outcome)

Q_stat_OV_DHEAS2<-mr_rucker(OV_DHEAS_2)
MR_PRESSO_OV_DHEAS2<-run_mr_presso(dat = OV_DHEAS_2)
MR_egger_OV_DHEAS2<-mr_egger_regression(b_exp = OV_DHEAS_2$beta.exposure,b_out = OV_DHEAS_2$beta.outcome,se_exp = OV_DHEAS_2$se.exposure,se_out = OV_DHEAS_2$se.outcome)

Q_stat_OV_ALDO2<-mr_rucker(OV_ALDO_2)
MR_PRESSO_OV_ALDO2<-run_mr_presso(dat = OV_ALDO_2)
MR_egger_OV_ALDO2<-mr_egger_regression(b_exp = OV_ALDO_2$beta.exposure,b_out = OV_ALDO_2$beta.outcome,se_exp = OV_ALDO_2$se.exposure,se_out = OV_ALDO_2$se.outcome)


##########################################################################Endometrial cancer#################################################################################
###Change this code for eight hormones:
UCEC_ANDRO<-harmonise_data(ANDRO_new,UCEC_GWAS,action = 2)
UCEC_TESTO<-harmonise_data(TESTO_new,UCEC_GWAS,action = 2)
UCEC_ESTRO<-harmonise_data(ESTRO_new,UCEC_GWAS,action = 2)
UCEC_PROG<-harmonise_data(PROG_new,UCEC_GWAS,action = 2)
UCEC_ALDO<-harmonise_data(ALDO_new,UCEC_GWAS,action = 2)
UCEC_CORT<-harmonise_data(CORT_new,UCEC_GWAS,action = 2)
UCEC_OHP17<-harmonise_data(OHP17_new,UCEC_GWAS,action = 2)
UCEC_DHEAS<-harmonise_data(DHEAS_new,UCEC_GWAS,action = 2)



UCEC_ANDRO_MR_Results<-mr(UCEC_ANDRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_TESTO_MR_Results<-mr(UCEC_TESTO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_ESTRO_MR_Results<-mr(UCEC_ESTRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_PROG_MR_Results<-mr(UCEC_PROG,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_CORT_MR_Results<-mr(UCEC_CORT,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_OHP17_MR_Results<-mr(UCEC_OHP17,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_DHEAS_MR_Results<-mr(UCEC_DHEAS,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_ALDO_MR_Results<-mr(UCEC_ALDO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))

Q_stat_UCEC_ANDRO<-mr_rucker(UCEC_ANDRO)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_ANDRO)
MR_egger_UCEC_ANDRO<-mr_egger_regression(b_exp = UCEC_ANDRO$beta.exposure,b_out = UCEC_ANDRO$beta.outcome,se_exp = UCEC_ANDRO$se.exposure,se_out = UCEC_ANDRO$se.outcome)

Q_stat_UCEC_TESTO<-mr_rucker(UCEC_TESTO)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_TESTO)
MR_egger_UCEC_TESTO<-mr_egger_regression(b_exp = UCEC_TESTO$beta.exposure,b_out = UCEC_TESTO$beta.outcome,se_exp = UCEC_TESTO$se.exposure,se_out = UCEC_TESTO$se.outcome)

Q_stat_UCEC_ESTRO<-mr_rucker(UCEC_ESTRO)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_ESTRO)
MR_egger_UCEC_ESTRO<-mr_egger_regression(b_exp = UCEC_ESTRO$beta.exposure,b_out = UCEC_ESTRO$beta.outcome,se_exp = UCEC_ESTRO$se.exposure,se_out = UCEC_ESTRO$se.outcome)

Q_stat_UCEC_PROG<-mr_rucker(UCEC_PROG)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_PROG)
MR_egger_UCEC_PROG<-mr_egger_regression(b_exp = UCEC_PROG$beta.exposure,b_out = UCEC_PROG$beta.outcome,se_exp = UCEC_PROG$se.exposure,se_out = UCEC_PROG$se.outcome)

Q_stat_UCEC_CORT<-mr_rucker(UCEC_CORT)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_CORT)
MR_egger_UCEC_CORT<-mr_egger_regression(b_exp = UCEC_CORT$beta.exposure,b_out = UCEC_CORT$beta.outcome,se_exp = UCEC_CORT$se.exposure,se_out = UCEC_CORT$se.outcome)

Q_stat_UCEC_OHP17<-mr_rucker(UCEC_OHP17)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_OHP17)
MR_egger_UCEC_OHP17<-mr_egger_regression(b_exp = UCEC_OHP17$beta.exposure,b_out = UCEC_OHP17$beta.outcome,se_exp = UCEC_OHP17$se.exposure,se_out = UCEC_OHP17$se.outcome)

Q_stat_UCEC_DHEAS<-mr_rucker(UCEC_DHEAS)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_DHEAS)
MR_egger_UCEC_DHEAS<-mr_egger_regression(b_exp = UCEC_DHEAS$beta.exposure,b_out = UCEC_DHEAS$beta.outcome,se_exp = UCEC_DHEAS$se.exposure,se_out = UCEC_DHEAS$se.outcome)

Q_stat_UCEC_ALDO<-mr_rucker(UCEC_ALDO)
MR_PRESSO_UCEC_AQNDRO<-run_mr_presso(dat = UCEC_ALDO)
MR_egger_UCEC_ALDO<-mr_egger_regression(b_exp = UCEC_ALDO$beta.exposure,b_out = UCEC_ALDO$beta.outcome,se_exp = UCEC_ALDO$se.exposure,se_out = UCEC_ALDO$se.outcome)



write.table(UCEC_ANDRO_MR_Results,file="UCEC_ANDRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(UCEC_ALDO_MR_Results,file="UCEC_ALDO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(UCEC_ESTRO_MR_Results,file="UCEC_ESTRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(UCEC_PROG_MR_Results,file="UCEC_PROG_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(UCEC_TESTO_MR_Results,file="UCEC_TESTO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(UCEC_CORT_MR_Results,file="UCEC_CORT_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(UCEC_OHP17_MR_Results,file="UCEC_OHP17_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(UCEC_DHEAS_MR_Results,file="UCEC_DHEAS_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)

#for bidirectional analysis
UCEC_GWAS_2<-UCEC_GWAS
names(UCEC_GWAS_2)<-gsub(x = names(UCEC_GWAS_2), pattern = "outcome", replacement = "exposure")
UCEC_GWAS_2<-UCEC_GWAS
names(UCEC_GWAS_2)<-gsub(x = names(UCEC_GWAS_2), pattern = "outcome", replacement = "exposure")
UCEC_GWAS_2<-UCEC_GWAS
names(UCEC_GWAS_2)<-gsub(x = names(UCEC_GWAS_2), pattern = "outcome", replacement = "exposure")
UCEC_GWAS_2<-UCEC_GWAS
names(UCEC_GWAS_2)<-gsub(x = names(UCEC_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")

ESTRO_new_2<-ESTRO_new
names(ESTRO_new_2)<-gsub(x = names(ESTRO_new_2), pattern = "exposure", replacement = "outcome")
PROG_new_2<-PROG_new
names(PROG_new_2)<-gsub(x = names(PROG_new_2), pattern = "exposure", replacement = "outcome")
TESTO_new_2<-TESTO_new
names(TESTO_new_2)<-gsub(x = names(TESTO_new_2), pattern = "exposure", replacement = "outcome")
ANDRO_new_2<-ANDRO_new
names(ANDRO_new_2)<-gsub(x = names(ANDRO_new_2), pattern = "exposure", replacement = "outcome")
ALDO_new_2<-ALDO_new
names(ALDO_new_2)<-gsub(x = names(ALDO_new_2), pattern = "exposure", replacement = "outcome")
OHP17_new_2<-OHP17_new
names(OHP17_new_2)<-gsub(x = names(OHP17_new_2), pattern = "exposure", replacement = "outcome")
DHEAS_new_2<-DHEAS_new
names(DHEAS_new_2)<-gsub(x = names(DHEAS_new_2), pattern = "exposure", replacement = "outcome")
CORT_new_2<-CORT_new
names(CORT_new_2)<-gsub(x = names(CORT_new_2), pattern = "exposure", replacement = "outcome")

UCEC_ANDRO_2<-harmonise_data(UCEC_GWAS_2,ANDRO_new_2,action = 2)
UCEC_TESTO_2<-harmonise_data(TESTO_new,UCEC_GWAS,action = 2)
UCEC_ESTRO_2<-harmonise_data(ESTRO_new,UCEC_GWAS,action = 2)
UCEC_PROG_2<-harmonise_data(PROG_new,UCEC_GWAS,action = 2)
UCEC_ALDO_2<-harmonise_data(ALDO_new,UCEC_GWAS,action = 2)
UCEC_CORT_2<-harmonise_data(CORT_new,UCEC_GWAS,action = 2)
UCEC_OHP17_2<-harmonise_data(OHP17_new,UCEC_GWAS,action = 2)
UCEC_DHEAS_2<-harmonise_data(DHEAS_new,UCEC_GWAS,action = 2)

UCEC_ANDRO_MR_Results2<-mr(UCEC_ANDRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_TESTO_MR_Results2<-mr(UCEC_TESTO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_ESTRO_MR_Results2<-mr(UCEC_ESTRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_PROG_MR_Results2<-mr(UCEC_PROG_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_CORT_MR_Results2<-mr(UCEC_CORT_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_OHP17_MR_Results2<-mr(UCEC_OHP17_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_DHEAS_MR_Results2<-mr(UCEC_DHEAS_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
UCEC_ALDO_MR_Results2<-mr(UCEC_ALDO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))


Q_stat_UCEC_ANDRO2<-mr_rucker(UCEC_ANDRO_2)
MR_PRESSO_UCEC_ANDRO2<-run_mr_presso(dat = UCEC_ANDRO_2)
MR_egger_UCEC_ANDRO2<-mr_egger_regression(b_exp = UCEC_ANDRO_2$beta.exposure,b_out = UCEC_ANDRO_2$beta.outcome,se_exp = UCEC_ANDRO_2$se.exposure,se_out = UCEC_ANDRO_2$se.outcome)

Q_stat_UCEC_TESTO2<-mr_rucker(UCEC_TESTO_2)
MR_PRESSO_UCEC_TESTO2<-run_mr_presso(dat = UCEC_TESTO_2)
MR_egger_UCEC_TESTO2<-mr_egger_regression(b_exp = UCEC_TESTO_2$beta.exposure,b_out = UCEC_TESTO_2$beta.outcome,se_exp = UCEC_TESTO_2$se.exposure,se_out = UCEC_TESTO_2$se.outcome)

Q_stat_UCEC_ESTRO2<-mr_rucker(UCEC_ESTRO_2)
MR_PRESSO_UCEC_ESTRO2<-run_mr_presso(dat = UCEC_ESTRO_2)
MR_egger_UCEC_ESTRO2<-mr_egger_regression(b_exp = UCEC_ESTRO_2$beta.exposure,b_out = UCEC_ESTRO_2$beta.outcome,se_exp = UCEC_ESTRO_2$se.exposure,se_out = UCEC_ESTRO_2$se.outcome)

Q_stat_UCEC_PROG2<-mr_rucker(UCEC_PROG_2)
MR_PRESSO_UCEC_PROG2<-run_mr_presso(dat = UCEC_PROG_2)
MR_egger_UCEC_PROG2<-mr_egger_regression(b_exp = UCEC_PROG_2$beta.exposure,b_out = UCEC_PROG_2$beta.outcome,se_exp = UCEC_PROG_2$se.exposure,se_out = UCEC_PROG_2$se.outcome)

Q_stat_UCEC_CORT2<-mr_rucker(UCEC_CORT_2)
MR_PRESSO_UCEC_CORT2<-run_mr_presso(dat = UCEC_CORT_2)
MR_egger_UCEC_CORT2<-mr_egger_regression(b_exp = UCEC_CORT_2$beta.exposure,b_out = UCEC_CORT_2$beta.outcome,se_exp = UCEC_CORT_2$se.exposure,se_out = UCEC_CORT_2$se.outcome)

Q_stat_UCEC_OHP172<-mr_rucker(UCEC_OHP17_2)
MR_PRESSO_UCEC_OHP172<-run_mr_presso(dat = UCEC_OHP17_2)
MR_egger_UCEC_OHP17<-mr_egger_regression(b_exp = UCEC_OHP17_2$beta.exposure,b_out = UCEC_OHP17_2$beta.outcome,se_exp = UCEC_OHP17_2$se.exposure,se_out = UCEC_OHP17_2$se.outcome)

Q_stat_UCEC_DHEAS2<-mr_rucker(UCEC_DHEAS_2)
MR_PRESSO_UCEC_DHEAS2<-run_mr_presso(dat = UCEC_DHEAS_2)
MR_egger_UCEC_DHEAS2<-mr_egger_regression(b_exp = UCEC_DHEAS_2$beta.exposure,b_out = UCEC_DHEAS_2$beta.outcome,se_exp = UCEC_DHEAS_2$se.exposure,se_out = UCEC_DHEAS_2$se.outcome)

Q_stat_UCEC_ALDO2<-mr_rucker(UCEC_ALDO_2)
MR_PRESSO_UCEC_ALDO2<-run_mr_presso(dat = UCEC_ALDO_2)
MR_egger_UCEC_ALDO2<-mr_egger_regression(b_exp = UCEC_ALDO_2$beta.exposure,b_out = UCEC_ALDO_2$beta.outcome,se_exp = UCEC_ALDO_2$se.exposure,se_out = UCEC_ALDO_2$se.outcome)

#for bidirectional analysis
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")
COLCA_GWAS_2<-COLCA_GWAS
names(COLCA_GWAS_2)<-gsub(x = names(COLCA_GWAS_2), pattern = "outcome", replacement = "exposure")

ESTRO_new_2<-ESTRO_new
names(ESTRO_new_2)<-gsub(x = names(ESTRO_new_2), pattern = "exposure", replacement = "outcome")
PROG_new_2<-PROG_new
names(PROG_new_2)<-gsub(x = names(PROG_new_2), pattern = "exposure", replacement = "outcome")
TESTO_new_2<-TESTO_new
names(TESTO_new_2)<-gsub(x = names(TESTO_new_2), pattern = "exposure", replacement = "outcome")
ANDRO_new_2<-ANDRO_new
names(ANDRO_new_2)<-gsub(x = names(ANDRO_new_2), pattern = "exposure", replacement = "outcome")
ALDO_new_2<-ALDO_new
names(ALDO_new_2)<-gsub(x = names(ALDO_new_2), pattern = "exposure", replacement = "outcome")
OHP17_new_2<-OHP17_new
names(OHP17_new_2)<-gsub(x = names(OHP17_new_2), pattern = "exposure", replacement = "outcome")
DHEAS_new_2<-DHEAS_new
names(DHEAS_new_2)<-gsub(x = names(DHEAS_new_2), pattern = "exposure", replacement = "outcome")
CORT_new_2<-CORT_new
names(CORT_new_2)<-gsub(x = names(CORT_new_2), pattern = "exposure", replacement = "outcome")

COLCA_ANDRO_2<-harmonise_data(COLCA_GWAS_2,ANDRO_new_2,action = 2)
COLCA_TESTO_2<-harmonise_data(TESTO_new,COLCA_GWAS,action = 2)
COLCA_ESTRO_2<-harmonise_data(ESTRO_new,COLCA_GWAS,action = 2)
COLCA_PROG_2<-harmonise_data(PROG_new,COLCA_GWAS,action = 2)
COLCA_ALDO_2<-harmonise_data(ALDO_new,COLCA_GWAS,action = 2)
COLCA_CORT_2<-harmonise_data(CORT_new,COLCA_GWAS,action = 2)
COLCA_OHP17_2<-harmonise_data(OHP17_new,COLCA_GWAS,action = 2)
COLCA_DHEAS_2<-harmonise_data(DHEAS_new,COLCA_GWAS,action = 2)

COLCA_ANDRO_MR_Results2<-mr(COLCA_ANDRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_TESTO_MR_Results2<-mr(COLCA_TESTO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_ESTRO_MR_Results2<-mr(COLCA_ESTRO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_PROG_MR_Results2<-mr(COLCA_PROG_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_CORT_MR_Results2<-mr(COLCA_CORT_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_OHP17_MR_Results2<-mr(COLCA_OHP17_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_DHEAS_MR_Results2<-mr(COLCA_DHEAS_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_ALDO_MR_Results2<-mr(COLCA_ALDO_2,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))


Q_stat_COLCA_ANDRO2<-mr_rucker(COLCA_ANDRO_2)
MR_PRESSO_COLCA_ANDRO2<-run_mr_presso(dat = COLCA_ANDRO_2)
MR_egger_COLCA_ANDRO2<-mr_egger_regression(b_exp = COLCA_ANDRO_2$beta.exposure,b_out = COLCA_ANDRO_2$beta.outcome,se_exp = COLCA_ANDRO_2$se.exposure,se_out = COLCA_ANDRO_2$se.outcome)

Q_stat_COLCA_TESTO2<-mr_rucker(COLCA_TESTO_2)
MR_PRESSO_COLCA_TESTO2<-run_mr_presso(dat = COLCA_TESTO_2)
MR_egger_COLCA_TESTO2<-mr_egger_regression(b_exp = COLCA_TESTO_2$beta.exposure,b_out = COLCA_TESTO_2$beta.outcome,se_exp = COLCA_TESTO_2$se.exposure,se_out = COLCA_TESTO_2$se.outcome)

Q_stat_COLCA_ESTRO2<-mr_rucker(COLCA_ESTRO_2)
MR_PRESSO_COLCA_ESTRO2<-run_mr_presso(dat = COLCA_ESTRO_2)
MR_egger_COLCA_ESTRO2<-mr_egger_regression(b_exp = COLCA_ESTRO_2$beta.exposure,b_out = COLCA_ESTRO_2$beta.outcome,se_exp = COLCA_ESTRO_2$se.exposure,se_out = COLCA_ESTRO_2$se.outcome)

Q_stat_COLCA_PROG2<-mr_rucker(COLCA_PROG_2)
MR_PRESSO_COLCA_PROG2<-run_mr_presso(dat = COLCA_PROG_2)
MR_egger_COLCA_PROG2<-mr_egger_regression(b_exp = COLCA_PROG_2$beta.exposure,b_out = COLCA_PROG_2$beta.outcome,se_exp = COLCA_PROG_2$se.exposure,se_out = COLCA_PROG_2$se.outcome)

Q_stat_COLCA_CORT2<-mr_rucker(COLCA_CORT_2)
MR_PRESSO_COLCA_CORT2<-run_mr_presso(dat = COLCA_CORT_2)
MR_egger_COLCA_CORT2<-mr_egger_regression(b_exp = COLCA_CORT_2$beta.exposure,b_out = COLCA_CORT_2$beta.outcome,se_exp = COLCA_CORT_2$se.exposure,se_out = COLCA_CORT_2$se.outcome)

Q_stat_COLCA_OHP172<-mr_rucker(COLCA_OHP17_2)
MR_PRESSO_COLCA_OHP172<-run_mr_presso(dat = COLCA_OHP17_2)
MR_egger_COLCA_OHP17<-mr_egger_regression(b_exp = COLCA_OHP17_2$beta.exposure,b_out = COLCA_OHP17_2$beta.outcome,se_exp = COLCA_OHP17_2$se.exposure,se_out = COLCA_OHP17_2$se.outcome)

Q_stat_COLCA_DHEAS2<-mr_rucker(COLCA_DHEAS_2)
MR_PRESSO_COLCA_DHEAS2<-run_mr_presso(dat = COLCA_DHEAS_2)
MR_egger_COLCA_DHEAS2<-mr_egger_regression(b_exp = COLCA_DHEAS_2$beta.exposure,b_out = COLCA_DHEAS_2$beta.outcome,se_exp = COLCA_DHEAS_2$se.exposure,se_out = COLCA_DHEAS_2$se.outcome)

Q_stat_COLCA_ALDO2<-mr_rucker(COLCA_ALDO_2)
MR_PRESSO_COLCA_ALDO2<-run_mr_presso(dat = COLCA_ALDO_2)
MR_egger_COLCA_ALDO2<-mr_egger_regression(b_exp = COLCA_ALDO_2$beta.exposure,b_out = COLCA_ALDO_2$beta.outcome,se_exp = COLCA_ALDO_2$se.exposure,se_out = COLCA_ALDO_2$se.outcome)



####################################################################Colorectal cancer#################################################################################
###Change this code for eight hormones:
COLCA_ANDRO<-harmonise_data(ANDRO_new,COLCA_GWAS,action = 2)
COLCA_TESTO<-harmonise_data(TESTO_new,COLCA_GWAS,action = 2)
COLCA_ESTRO<-harmonise_data(ESTRO_new,COLCA_GWAS,action = 2)
COLCA_PROG<-harmonise_data(PROG_new,COLCA_GWAS,action = 2)
COLCA_ALDO<-harmonise_data(ALDO_new,COLCA_GWAS,action = 2)
COLCA_CORT<-harmonise_data(CORT_new,COLCA_GWAS,action = 2)
COLCA_OHP17<-harmonise_data(OHP17_new,COLCA_GWAS,action = 2)
COLCA_DHEAS<-harmonise_data(DHEAS_new,COLCA_GWAS,action = 2)



COLCA_ANDRO_MR_Results<-mr(COLCA_ANDRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_TESTO_MR_Results<-mr(COLCA_TESTO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_ESTRO_MR_Results<-mr(COLCA_ESTRO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_PROG_MR_Results<-mr(COLCA_PROG,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_CORT_MR_Results<-mr(COLCA_CORT,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_OHP17_MR_Results<-mr(COLCA_OHP17,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_DHEAS_MR_Results<-mr(COLCA_DHEAS,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
COLCA_ALDO_MR_Results<-mr(COLCA_ALDO,parameters = default_parameters(),method_list =c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))

Q_stat_COLCA_ANDRO<-mr_rucker(COLCA_ANDRO)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_ANDRO)
MR_egger_COLCA_ANDRO<-mr_egger_regression(b_exp = COLCA_ANDRO$beta.exposure,b_out = COLCA_ANDRO$beta.outcome,se_exp = COLCA_ANDRO$se.exposure,se_out = COLCA_ANDRO$se.outcome)

Q_stat_COLCA_TESTO<-mr_rucker(COLCA_TESTO)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_TESTO)
MR_egger_COLCA_TESTO<-mr_egger_regression(b_exp = COLCA_TESTO$beta.exposure,b_out = COLCA_TESTO$beta.outcome,se_exp = COLCA_TESTO$se.exposure,se_out = COLCA_TESTO$se.outcome)

Q_stat_COLCA_ESTRO<-mr_rucker(COLCA_ESTRO)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_ESTRO)
MR_egger_COLCA_ESTRO<-mr_egger_regression(b_exp = COLCA_ESTRO$beta.exposure,b_out = COLCA_ESTRO$beta.outcome,se_exp = COLCA_ESTRO$se.exposure,se_out = COLCA_ESTRO$se.outcome)

Q_stat_COLCA_PROG<-mr_rucker(COLCA_PROG)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_PROG)
MR_egger_COLCA_PROG<-mr_egger_regression(b_exp = COLCA_PROG$beta.exposure,b_out = COLCA_PROG$beta.outcome,se_exp = COLCA_PROG$se.exposure,se_out = COLCA_PROG$se.outcome)

Q_stat_COLCA_CORT<-mr_rucker(COLCA_CORT)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_CORT)
MR_egger_COLCA_CORT<-mr_egger_regression(b_exp = COLCA_CORT$beta.exposure,b_out = COLCA_CORT$beta.outcome,se_exp = COLCA_CORT$se.exposure,se_out = COLCA_CORT$se.outcome)

Q_stat_COLCA_OHP17<-mr_rucker(COLCA_OHP17)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_OHP17)
MR_egger_COLCA_OHP17<-mr_egger_regression(b_exp = COLCA_OHP17$beta.exposure,b_out = COLCA_OHP17$beta.outcome,se_exp = COLCA_OHP17$se.exposure,se_out = COLCA_OHP17$se.outcome)

Q_stat_COLCA_DHEAS<-mr_rucker(COLCA_DHEAS)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_DHEAS)
MR_egger_COLCA_DHEAS<-mr_egger_regression(b_exp = COLCA_DHEAS$beta.exposure,b_out = COLCA_DHEAS$beta.outcome,se_exp = COLCA_DHEAS$se.exposure,se_out = COLCA_DHEAS$se.outcome)

Q_stat_COLCA_ALDO<-mr_rucker(COLCA_ALDO)
MR_PRESSO_COLCA_AQNDRO<-run_mr_presso(dat = COLCA_ALDO)
MR_egger_COLCA_ALDO<-mr_egger_regression(b_exp = COLCA_ALDO$beta.exposure,b_out = COLCA_ALDO$beta.outcome,se_exp = COLCA_ALDO$se.exposure,se_out = COLCA_ALDO$se.outcome)


write.table(COLCA_ANDRO_MR_Results,file="COLCA_ANDRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(COLCA_ALDO_MR_Results,file="COLCA_ALDO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(COLCA_ESTRO_MR_Results,file="COLCA_ESTRO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(COLCA_PROG_MR_Results,file="COLCA_PROG_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(COLCA_TESTO_MR_Results,file="COLCA_TESTO_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(COLCA_CORT_MR_Results,file="COLCA_CORT_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(COLCA_OHP17_MR_Results,file="COLCA_OHP17_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(COLCA_DHEAS_MR_Results,file="COLCA_DHEAS_MR_Results.txt",quote=FALSE,sep="\t",row.names=FALSE)


######

res_single <- mr_singlesnp(dat = PRAD_ANDRO)
p2 <- mr_forest_plot(res_single)

res_single <- mr_singlesnp(dat=UCEC_ANDRO)
p4 <- mr_funnel_plot(res_single)
p4[[1]]