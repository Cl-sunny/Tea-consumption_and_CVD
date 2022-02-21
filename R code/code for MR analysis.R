rm(list=ls())
library(TwoSampleMR)
library(MRPRESSO)
library(meta)
exposure_dat <- read_exposure_data(
  filename = 'SNP-tea consumption associations.csv',
  sep = ',',
  snp_col = 'SNPs',
  beta_col = 'beta',
  se_col = 'se',
  effect_allele_col = 'ALT',
  phenotype_col = 'tea consumption',
  other_allele_col = 'REF',
  eaf_col = 'minor_AF',
  pval_col = 'pval'
)
exposure_dat$exposure <- 'tea consumption'
outcome_dat <- read.csv("SNP-CVD associations.csv")
dat <- harmonise_data(exposure_dat,outcome_dat, action = 1)
#####P-value cut-off 5E-8###
dat1=dat[dat$pval.exposure<5e-8,]
mr_results1 <- mr(dat1,method_list=c("mr_ivw_mre",
                                     "mr_weighted_median",
                                     "mr_egger_regression"))
het1=mr_heterogeneity(dat1)
ple1=mr_pleiotropy_test(dat1)
OR_mr_results1=generate_odds_ratios(mr_results1)

###MR-PRESSO ANALYSIS###
mr_presso1=run_mr_presso(dat1, NbDistribution = 2500, SignifThreshold = 0.05)

######META-ANALYSIS###
data_meta=mr_results1[mr_results1$method=="Inverse variance weighted (multiplicative random effects)",]
meta_CAD=data_meta[grep("CAD",data_meta$outcome),]
m_CAD <-metagen(meta_CAD$b,meta_CAD$se,studlab=meta_CAD$outcome, sm="BETA",backtransf=TRUE)
m_CAD

meta_MI=data_meta[grep("MI",data_meta$outcome),]
m_MI <-metagen(meta_MI$b,meta_MI$se,studlab=meta_MI$outcome, sm="BETA",backtransf=TRUE)
m_MI

meta_AF=data_meta[grep("AF",data_meta$outcome),]
m_AF <-metagen(meta_AF$b,meta_AF$se,studlab=meta_AF$outcome, sm="BETA",backtransf=TRUE)
m_AF

meta_HF=data_meta[grep("HF",data_meta$outcome),]
m_HF <-metagen(meta_HF$b,meta_HF$se,studlab=meta_HF$outcome, sm="BETA",backtransf=TRUE)
m_HF


#####P-value cut-off 5E-7###
dat2=dat[dat$pval.exposure<5e-7,]
mr_results2 <- mr(dat2,method_list=c("mr_ivw_mre",
                                     "mr_weighted_median",
                                     "mr_egger_regression"))
het2=mr_heterogeneity(dat2)
ple2=mr_pleiotropy_test(dat2)
OR_mr_results2=generate_odds_ratios(mr_results2)

###MR-PRESSO ANALYSIS###
mr_presso2=run_mr_presso(dat2, NbDistribution = 2500, SignifThreshold = 0.05)

######META-ANALYSIS###
data_meta=mr_results2[mr_results2$method=="Inverse variance weighted (multiplicative random effects)",]
meta_CAD=data_meta[grep("CAD",data_meta$outcome),]
m_CAD <-metagen(meta_CAD$b,meta_CAD$se,studlab=meta_CAD$outcome, sm="BETA",backtransf=TRUE)
m_CAD

meta_MI=data_meta[grep("MI",data_meta$outcome),]
m_MI <-metagen(meta_MI$b,meta_MI$se,studlab=meta_MI$outcome, sm="BETA",backtransf=TRUE)
m_MI

meta_AF=data_meta[grep("AF",data_meta$outcome),]
m_AF <-metagen(meta_AF$b,meta_AF$se,studlab=meta_AF$outcome, sm="BETA",backtransf=TRUE)
m_AF

meta_HF=data_meta[grep("HF",data_meta$outcome),]
m_HF <-metagen(meta_HF$b,meta_HF$se,studlab=meta_HF$outcome, sm="BETA",backtransf=TRUE)
m_HF
