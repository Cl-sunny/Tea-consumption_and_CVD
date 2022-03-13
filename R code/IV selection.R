rm(list=ls())
library(data.table) 
Allvariants <- fread('1488_raw.gwas.imputed_v3.both_sexes.tsv',header=T) 
head(Allvariants) 

#select SNPs with P-values <5e-7
selected_variant <- Allvariants[Allvariants$pval<5e-7,] 

#select SNPs with minor frequencies > 1%
selected_variant <- selected_variant[selected_variant$minor_AF>0.01,] 

#select SNPs autosomal biallelic SNPs
library(stringr)
selected_variant$CHR <- str_split(as.character(selected_variant$variant),':',simplify = T)[,1]
selected_variant$POS <- str_split(as.character(selected_variant$variant),':',simplify = T)[,2]
selected_variant$REF <- str_split(as.character(selected_variant$variant),':',simplify = T)[,3]
selected_variant$ALT <- str_split(as.character(selected_variant$variant),':',simplify = T)[,4]
selected_variant <- selected_variant[!selected_variant$CHR%in%c("X","23"),]
selected_variant <- selected_variant[nchar(selected_variant$REF)<2,]
selected_variant <- selected_variant[nchar(selected_variant$ALT)<2,]
selected_variant$CHR_POS <- paste(selected_variant$CHR,selected_variant$POS,sep=":")
selected_variant[duplicated(selected_variant$CHR_POS),] 
###there was a duplicated variant 15:75139426 and we keep the SNP with the smaller P-value
selected_variant <- selected_variant[order(selected_variant$pval),]
selected_variant <- selected_variant[!duplicated(selected_variant$CHR_POS),]

#write.csv(selected_variant,"selected_variant.csv")

###As the summary statistics for tea consumption from Neale lab only provided the information on chromosome and position, we have to find out rs number for each genetic variant.
selected_variant <- read.csv("information on selected variant.csv")

##perform clumping for SNPs with P-values<5e-8
library(ieugwasr)
mydata <- selected_variant[selected_variant$pval<5e-8,c("SNPs","pval")]
colnames(mydata) <- c("rsid","pval")
IV1 <- ld_clump(dat = mydata,clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99,
                 pop = "EUR", access_token = NULL,bfile = NULL, plink_bin = NULL
) 
data_IV1 <- selected_variant[selected_variant$SNPs%in%IV1$rsid,c("SNPs","CHR","POS",
                                                                 "ALT","REF","minor_AF",
                                                                 "beta","se","pval")]


##perform clumping for SNPs with P-values<5e-7
mydata2 <- selected_variant[selected_variant$pval<5e-7,c("SNPs","pval")]
colnames(mydata2) <- c("rsid","pval")
IV2 <- ld_clump(dat = mydata2,clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99,
                pop = "EUR", access_token = NULL,bfile = NULL, plink_bin = NULL
) 
data_IV2 <- selected_variant[selected_variant$SNPs%in%IV2$rsid,c("SNPs","CHR","POS",
                                                                 "ALT","REF","minor_AF",
                                                                 "beta","se","pval")]

