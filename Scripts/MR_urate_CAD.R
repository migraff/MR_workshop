#title: "MR_urate_CAD"

#tell R to use the packages you installed


###--- install Mendelian randomization package: "MendelianRandomization" ---###
###--- install Meta-Analysis Package: "metafor" - used here for forest plots---###
###--- install devtools to install from github, etc
###--- install MRbase items; "MRCIEU/TwoSampleMR" - used here for LD pruning functions ---###
###--- install MRPRESSO: MR method ---###
###--- install tidyverse: helps streamline coding --###
###--- install data.table: here to read in large gz files --###


###--- install Mendelian randomization package: "MendelianRandomization" ---###
#install.packages("MendelianRandomization", dependencies = TRUE)
library(MendelianRandomization)


###--- install Meta-Analysis Package: "metafor" - used here for forest plots---###
#install.packages("metafor", dependencies = TRUE)
library(metafor)

###--- install devtools to install from github, etc
#install.packages("devtools")
library(devtools)

###--- install MRbase items; "MRCIEU/TwoSampleMR" - used here for LD pruning functions ---###
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)


###--- install MRPRESSO: MR method ---###
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")

library(MRPRESSO)


##-- install tidyverse: helps streamline coding --###
#install.packages("tidyverse", dependencies = TRUE)

library(tidyverse)

###--- install data.table: to read in large gz files --###
#install.packages("data.table", dependencies = TRUE)

library(data.table)


#Set working directory

setwd("/proj/epi/CVDGeneNas/migraff/MR/Workshop/Nov2023/")

#1 i.  read in GWAS for exposure file, a GWAS of uric acid.  Original file is called "SOL_urate_summary.csv.gz".
	##For this workshop I've reduced to a smaller file of variants, those with p<0.01.
	##We will limit this file to variants that have pvalue<5x10-6

#original exposure file
#sol_urate<-read_csv(gzfile("Exposure_GWAS/SOL_urate_summary.csv.gz"),col_name=TRUE)

#exposure file for workshop
sol_urate<-read_csv(gzfile("Exposure_GWAS_SOL_urate.csv.gz"),col_name=TRUE)
sol_urate_lim<-sol_urate[sol_urate$PVAL<0.000005,] 
nrow(sol_urate_lim)


#1 ii. create a variant ID with chr, pos, and alleles with urate GWAS to merge with outcome

sol_urate_lim <- sol_urate_lim %>% 
	mutate(cpaid = str_c(chromosome,position,pmin(alleleA, alleleB),
		pmax(alleleA,alleleB), sep=":"))

#1 iii. rename the columns in the exposure file to be the same for every file

exposure<-data.frame(cpaid = sol_urate_lim$cpaid,
                  beta.exposure = sol_urate_lim$Beta,
                  se.exposure = sol_urate_lim$SE,
                  effect_allele.exposure = sol_urate_lim$alleleA,
                  other_allele.exposure =sol_urate_lim$alleleB,
                  eaf.exposure = sol_urate_lim$AF,
                  pval.exposure = sol_urate_lim$PVAL)
exposure[1,]


#2 i. read downloadeded Cardiogram file - GWAS of Coronary Artery Disease (CAD) 
#Original file is called "cad.add.160614.website.txt.gz".
	##For this workshop I've reduced to a smaller file of variants that overlap with variants selected from exposure file above.
	##We will limit this file to variants that have pvalue<5x10-6

 
#Original GWAS file
#cad<-fread("Outcome_GWAS/cad.add.160614.website.txt.gz")

#Outcome GWAS file for workshop
cad<-read.table("Outcome_GWAS_CAD_160614.txt", header=T, sep="\t")
nrow(cad)

 

#2 ii. create a variant ID with chr, pos, and alleles for merging with exposure
 
cad <- cad %>% mutate(cpaid = str_c(chr,bp_hg19,pmin(effect_allele,noneffect_allele),
			pmax(effect_allele,noneffect_allele), sep=":"))

 

#2 iii. select the genetic instrument SNPs from the outcome file 
 
sl<-cad[cad$cpaid %in% exposure$cpaid, ] # get variants in CAD GWAS
nrow(sl)

 

#2 iv. rename the columns in the outcome file to be the same for every file
 
outcome<-data.frame(cpaid = sl$cpaid,
                  beta.outcome = sl$beta,
                  se.outcome = sl$se_dgc,
                  effect_allele.outcome = sl$effect_allele,
                  other_allele.outcome =sl$noneffect_allele,
                  eaf.outcome = sl$effect_allele_freq,
                  pval.outcome = sl$p_dgc,
		  SNP=sl$markername)
outcome[1,]
 

#3 Merge SNP to exposure and SNP to outcome files
 
exp_out<-merge(outcome,exposure,by="cpaid",all.x=TRUE) #merge SNP(cpaid) to exposure and SNP (cpaid) to outcome file
nrow(exp_out)
exp_out[3,]
 

#4. use MRbase to find indepnednt SNPs with r2<0.05
#then detach package because it precludes the IVW to work in the MendelianRandomization package
 
all<-clump_data(exp_out,clump_r2=0.05,pop="AMR")

nrow(all)
detach(package:TwoSampleMR)
 

#5 calculate F-statistic for each SNP, then look at the mean
 
all$f1<-(all$beta.exposure*all$beta.exposure)/(all$se.exposure*all$se.exposure)
all$f1
mean(all$f1)
 

#6 Align the SNPs on the same effect allele for exposure and outcome
# by changing sign of beta.exposure if effect alleles do not match
 
all$effect_allele.outcome <- as.factor(all$effect_allele.outcome)
all$effect_allele.exposure <- as.factor(all$effect_allele.exposure)
lev2 <- unique( c( levels(all$effect_allele.outcome), levels(all$effect_allele.exposure) ) )
all$effect_allele.outcome <- factor(all$effect_allele.outcome, levels=lev2)
all$effect_allele.exposure <- factor(all$effect_allele.exposure, levels=lev2)
all$effect_allele.exposure<-gsub(" ", "",all$effect_allele.exposure, fixed = TRUE)
all$beta.exposure[all$effect_allele.exposure!=all$effect_allele.outcome]<-all$beta.exposure[all$effect_allele.exposure!=all$effect_allele.outcome] * -1
all[3,]
 

#7 .get forest plot with fixed effects, same as Mendelianrandomization IVW with fixed effects 
 
x<-all$beta.exposure # beta for SNP to exposure
sigmax<-all$se.exposure # its standard errors
y<-all$beta.outcome # beta for SNP to outcome
sigmay<-all$se.outcome # its standard errors
all$Wald<-y/x #Wald estimate
all$Waldvar<-(sigmay^2/x^2) # using Burgess's method
#all$lab<-paste(all$SNP, sep=" ")
dmres<-rma.uni(yi=all$Wald, vi=all$Waldvar, slab=all$SNP, method="FE")
dmres

par(mfrow=c(2,2)) 
par(mar=c(4,2,0,2)+.1)
png(file="forestplot_sol_urate_CAD.png")
forest(dmres, atransf=exp,xlab=" ", mlab="Coronary Artery disease (OR)", at=log(c(.5, 1,2)),xlim=c(-1.7,1.3),cex=.8)
title("MR uric acid to CAD")
dev.off()
 



#8i. save the data for future use (RESCUE FILE)

 
write.csv(all,"urate_CAD.csv")

 

#8 ii uncomment and read in file of SNPs to exposure and outcome for lead if you have problems with the above
##this is the main data needed to run the MR analyses below

 
#all<-read.csv("../urate_CAD.csv",header=TRUE)
#nrow(all)
 

#9 get estimates for ihd using Mendelian randomization package with fixed effects
#nb exponentiate to get OR and 95% CIs matching the forest plot

 
MRInputObject <- mr_input(all$beta.exposure, all$se.exposure, all$beta.outcome, all$se.outcome)
mr_ivw(MRInputObject,model="fixed")

 

#10 get MR estimates for ihd using random effects
 
mr_ivw(MRInputObject)
 

#11 get weighted median and MR-Egger estimates
 
mr_median(MRInputObject)
mr_egger(MRInputObject)
 

#12 run MR-Presso; this method helps identify snps that are outliers as part of the IV 
# 12 i first we will run the simulated toy dataset to show
 

data(SummaryStats)

# Run MR-PRESSO global method
mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = SummaryStats, NbDistribution = 1000,  SignifThreshold = 0.05)


#identify the SNPs MR-Pressos said were outliers
SummaryStats[7,]$SNP
 
# 12 ii now load our data with uric acid and CAD to check for outliers

 
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = all, NbDistribution = 1000,  SignifThreshold = 0.05)

 

#lets have a look at the original forest plot again
#4 get forest plot with fixed effects, same as Mendelianrandomization IVW with fixed effects 
 
x<-all$beta.exposure # beta for SNP to exposure
sigmax<-all$se.exposure # its standard errors
y<-all$beta.outcome # beta for SNP to outcome
sigmay<-all$se.outcome # its standard errors
all$Wald<-y/x #Wald estimate
all$Waldvar<-(sigmay^2/x^2) # using Burgess's method
all$lab<-paste(all$SNP, all$gene, sep=" ")
#all<-all[all$SNP!="rs550057",] #uncomment this row to exclude rs550057
dmres<-rma.uni(yi=all$Wald, vi=all$Waldvar, slab=all$lab, method="FE")
dmres
dmres<-rma.uni(yi=all$Wald, vi=all$Waldvar, slab=all$lab, method="REML")
dmres
forest(dmres, atransf=exp,xlab=" ", mlab="Coronary Artery disease (OR)", at=log(c(.5, 1,2)),xlim=c(-1.7,1.3),cex=.8)

