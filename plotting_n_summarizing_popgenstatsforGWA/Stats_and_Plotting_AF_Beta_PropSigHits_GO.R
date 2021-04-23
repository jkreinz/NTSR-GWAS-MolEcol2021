
library(data.table)
library(dplyr)
library(ggplot2)
library(PNWColors)

########################
#AF stats and plotting
########################


library(data.table)
library(dplyr)
library(ggplot2)

#read in and organize
res<-fread("~/nonbinary_twocovs_n155.assoc.txt",header=T) #gemma observed GWA results
res$fdr_p<-p.adjust(res$p_lrt,method="fdr")
names(res)[c(1,3)]<-c("CHROM","POS")

#read in and organize
permuted_sigSNPS<-read.table("twocovs_FDR05_sigsnps.txt",header = F) #gemma within-pop permuted GWA results
names(permuted_sigSNPS)[c(1,2,4,8)]<-c("itt","CHROM","POS","af")
permuted_sigSNPS<-permuted_sigSNPS[,c(2,4,8)]

#read in and organize
permuted_region_sigSNPS<-read.table("twocovs_regional_FDR05_sigsnps.txt",header = F) #gemma within-region permuted GWA results
names(permuted_region_sigSNPS)[c(1,2,4,8)]<-c("itt","CHROM","POS","af")
permuted_region_sigSNPS<-permuted_region_sigSNPS[,c(2,4,8)]

#make class subsets based on signficance
TSRsig<-res[res$fdr_p<0.05,]
TSRsig<-TSRsig[,c(1,3,7)]

TSRsig_not<-anti_join(res,TSRsig,by=c("CHROM","POS"))
TSRsig_not<-TSRsig_not[,c(1,3,7)]

TSRsig_not$Class<-"Insig."
TSRsig$Class<-"Sig."
permuted_sigSNPS$Class<-"Pop Permuted"
permuted_region_sigSNPS$Class<-"Region Permuted"

#merge classes into long format data
allmerged <- do.call("rbind", list(TSRsig_not, permuted_region_sigSNPS, permuted_sigSNPS, TSRsig))

#establish pretty color scheme
bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[1],bay[3],bay[6],bay[8])

#relevel factors
allmerged$Class <- factor(allmerged$Class, levels=c("Insig.","Pop Permuted","Region Permuted","Sig."))

#plot long data grouped by class
ggplot(data=allmerged, aes(af)) +
  geom_density(alpha=.6, aes(y=..scaled..,color=Class,fill=Class)) +
  theme_classic() +
  scale_color_manual(values = three) +
  scale_fill_manual(values = three) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Frequency",x="Allele frequency")

#get AF stats - bootstrap sample permuted set
x=nrow(permuted_sigSNPS)
y=573

perm_resample<-list(); perm_mean<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-permuted_sigSNPS[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_median[[i]] <- median(na.omit(interm$af))
  perm_variance[[i]] <- var(na.omit(interm$af))
  
}

#make dataframe of bootstrapped results
perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_median))
names(perm_stats)<-"median_AF"
perm_stats$variance<-unlist(perm_variance)

quantile(perm_stats$median_AF,probs = c(0.05, 0.95)) #95% median CIs of pop-permuted GWA
median(TSRsig$af) #median of observed GWA



########################
#BETA stats and plotting
########################


library(data.table)
library(dplyr)
library(ggplot2)

#read in and organize
res<-fread("~/nonbinary_twocovs_n155.assoc.txt",header=T) #gemma observed GWA results
res$fdr_p<-p.adjust(res$p_lrt,method="fdr")
names(res)[c(1,3)]<-c("CHROM","POS")

#read in and organize
permuted_sigSNPS<-read.table("twocovs_FDR05_sigsnps.txt",header = F) #gemma within-pop permuted GWA results
names(permuted_sigSNPS)[c(1,2,4,8)]<-c("itt","CHROM","POS","af")
permuted_sigSNPS<-permuted_sigSNPS[,c(2,4,8)]

#read in and organize
permuted_region_sigSNPS<-read.table("twocovs_regional_FDR05_sigsnps.txt",header = F) #gemma within-region permuted GWA results
names(permuted_region_sigSNPS)[c(1,2,4,8)]<-c("itt","CHROM","POS","af")
permuted_region_sigSNPS<-permuted_region_sigSNPS[,c(2,4,8)]

#make class subsets based on signficance
TSRsig<-res[res$fdr_p<0.05,]
TSRsig<-TSRsig[,c(1,3,7)]

TSRsig_not<-anti_join(res,TSRsig,by=c("CHROM","POS"))
TSRsig_not<-TSRsig_not[,c(1,3,7)]

TSRsig_not$Class<-"Insig."
TSRsig$Class<-"Sig."
permuted_sigSNPS$Class<-"Pop Permuted"
permuted_region_sigSNPS$Class<-"Region Permuted"

#subset columns wrt to BETA
permuted_sigSNPS<-permuted_sigSNPS[,c(2,4,8,9,10)]
permuted_sigSNPS$Class<-"Pop Permuted"
permuted_regional_sigSNPS<-permuted_regional_sigSNPS[,c(2,4,8,9,10)]
permuted_regional_sigSNPS$Class<-"Region Permuted"

names(permuted_sigSNPS)<-names(TSRsig) #carry over names
names(permuted_regional_sigSNPS)<-names(TSRsig) #carry over names

both<-rbind(permuted_sigSNPS,TSRsig) #merge
both<-rbind(both,permuted_regional_sigSNPS)

#set up SEs for error bars
both$upper<-both$beta+both$se
both$lower<-both$beta-both$se
TSRsig$lower<-TSRsig$beta-TSRsig$se
TSRsig$upper<-TSRsig$beta+TSRsig$se
permuted_sigSNPS$lower<-permuted_sigSNPS$beta-permuted_sigSNPS$se
permuted_sigSNPS$upper<-permuted_sigSNPS$beta+permuted_sigSNPS$se
permuted_regional_sigSNPS$lower<-permuted_regional_sigSNPS$beta-permuted_regional_sigSNPS$se
permuted_regional_sigSNPS$upper<-permuted_regional_sigSNPS$beta+permuted_regional_sigSNPS$se


bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[3],bay[6],bay[8])

#relevel
both$Class <- factor(both$Class, levels=c("Permuted","Observed"))

#plot AF by Beta, for observed and permuted sets
ggplot() +
  geom_errorbar(data=permuted_sigSNPS, aes(x=af, ymax=upper,ymin=lower),alpha=0.05,color=bay[6] ) +
  geom_point(data=permuted_sigSNPS,aes(af,beta),alpha=.02,color=bay[6],cex=3) +
  geom_errorbar(data=permuted_regional_sigSNPS, aes(x=af, ymax=upper,ymin=lower),alpha=0.05,color=bay[3] ) +
  geom_point(data=permuted_regional_sigSNPS,aes(af,beta),alpha=.02,color=bay[3],cex=3) +
  geom_errorbar(data=TSRsig, aes(x=af, ymax=upper,ymin=lower),alpha=0.3,color=bay[8] ) +
  geom_point(data=TSRsig,aes(af,beta),alpha=.4,cex=3, color=bay[8]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=three) +
  scale_color_manual(values=three) +
  labs(y="Effect size (Beta)",x="Allele frequency") +
  ylim(-3,4)

#calculate bootstrap stats
x <- nrow(permuted_sigSNPS)
y <- 573

perm_resample<-list(); perm_median<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-permuted_sigSNPS[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  neg<-sum(interm$beta<0) #count number of sig negative effect sizes
  pos<-sum(interm$beta>=0)  #count number of sig positive effect sizes
  propneg_beta[[i]] <-neg/(pos+neg) #proportion of negative effects
  
}

#make dataframe out of bootstrapped results
perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, propneg_beta))
names(perm_stats)<-"propneg_beta"

quantile(perm_stats$propneg_beta,probs = c(0.05, 0.95)) #permuted 95% CIs of proportion of negative effects
neg<-sum(TSRsig$beta<0) #observed # of negative effects
pos<-sum(TSRsig$beta>=0)  #observed # of positive effects
neg/(pos+neg) #observed prop neg effects

#now calculate expected abs effect size from bootstrapping
x <- nrow(permuted_sigSNPS)
y <- 573

perm_resample<-list(); perm_median<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-permuted_sigSNPS[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_median[[i]] <- median(abs(interm$beta))
  
}

perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_median))
names(perm_stats)<-"median_beta_abs"
#perm_stats$variance<-unlist(perm_variance)

quantile(perm_stats$median_beta_abs,probs = c(0.05, 0.95)) #expected CIs
median(abs(TSRsig$beta)) #observed

#############################################
#Number of sig hits by gene and chromosome
#############################################


####################
#get number of significant hits by gene
####################

hitsbygene<-read.table("nonbinary_FDR05_observedsigsnps.gff",fill = T,na.strings = "") #GFF filtered to include genes w sig hits, two covariates
hitsbygene<-hitsbygene[,1:12]
hitsbygene<-hitsbygene[complete.cases(hitsbygene),]

library(tidyr)
hitsbygene<-hitsbygene %>% separate(V12,",", extra = "drop", fill = "right") #more cleaning

gwas_TSR<-hitsbygene %>% group_by(V1) %>% tally #summarise hits by gene name
gwas_TSR$V1 <- factor(gwas_TSR$V1, levels = c("Scaffold_1","Scaffold_2","Scaffold_3","Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"))

#and import no covariate genes
hitsbygene_nocov<-read.table("NTSR_vanilla_n155.gff",fill = T, na.strings = "")
hitsbygene_nocov<-hitsbygene_nocov[,1:13]
hitsbygene_nocov<-hitsbygene_nocov[complete.cases(hitsbygene_nocov),]
hitsbygene_nocov <- hitsbygene_nocov[hitsbygene_nocov$V10 == ".",]

gwas_nocov<-hitsbygene_nocov %>% group_by(V1) %>% tally
gwas_nocov$V1 <- factor(gwas_nocov$V1, levels = c("Scaffold_1","Scaffold_2","Scaffold_3","Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"))

#make new variable called model (referring to GWA model)
gwas_TSR$model<-"bothTSR_covs"
gwas_nocov$model<-"No_covariates"
gwas_TSR$prop<-gwas_TSR$n/(sum(gwas_TSR$n))
gwas_nocov$prop<-gwas_nocov$n/(sum(gwas_nocov$n))

gwas_both<-rbind(gwas_TSR,gwas_nocov)

#install.packages("PNWColors")
library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[4],bay[2])

library(tidyverse)
library(forcats)
four$model <- factor(four$model, levels = c("No_covariates","bothTSR_covs"))

#plot prportion of significant hits per scaffold by model
four %>% filter(V1 != "/") %>% mutate(V1 = str_remove_all(V1, "Scaffold_")) %>% 
  mutate(V1 = factor(V1, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))) %>%
  ggplot(aes(V1,prop, color=model)) +
  geom_point(cex=3,alpha=.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))  +
  labs(x="Scaffold",y="Proportion of Sig. SNPs") +
  scale_color_manual(values=c("grey60","grey30"))

####################
#PLOT GO ENRICHMENT
####################

go<-read.csv("Desktop/NTSR_enrichment_redo.csv",header=T)
library(ggplot2)
library(stringr)
#go$bonferonni[go$bonferonni == 0] <- ""
#go$bonferonni[go$bonferonni == 1] <- "*"

#justbon<-go[go$bonferonni==1,]

ggplot(data=go, aes(observed,reorder(str_to_title(GO.biological.process.complete),(-raw.P.value)), fill=as.numeric(class))) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c() +
  ylab("GO Biological Processes") +
  xlab("Number of Observed Genes") +
  theme_bw() +
  xlim(0,60) +
  geom_text(aes(label=bonferonni,hjust=-.2,vjust=.7)) +
  guides(fill=guide_legend(title="Heiarchcal \ncategory"))


