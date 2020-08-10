library(ggplot2)
library(tidyverse)
library(PNWColors)
library(data.table)

#read in observed, expected (permuted), and genomewide sets
ihs_obs<-read.table("allobservedhits_IHS.txt")
ihs_exp<-read.table("allpermutedhits_IHS.txt")
ihs_gw<-fread("allNOTobservedhits_IHS.txt")

names(ihs_obs)<-c("snp","ps","freq","iHS1","iHS0","unstand","left","right","ancestral_left","ancestral_right")
names(ihs_exp)<-c("snp","ps","freq","iHS1","iHS0","unstand","left","right","ancestral_left","ancestral_right")
names(ihs_gw)<-c("snp","ps","freq","iHS1","iHS0","unstand","left","right","ancestral_left","ancestral_right")

#swap iHH1 and iHH0 values such that all 1 alleles have positive effects on resistance
sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
hdr<-names(sex)
sex$FDR_corr<-p.adjust(sex$p_lrt,method = "fdr")
sig<-sex[sex$FDR_corr<0.05,]
sig_neg<-sig[sig$beta<0,]

neg_effects<-semi_join(ihs_obs,sig_neg,by=c("ps"))
pos_effects<-anti_join(ihs_obs,sig_neg,by=c("ps"))
neg_swapped<-neg_effects[,c(1,2,3,5,4,6,7,8,9,10)]
obs<-rbind(pos_effects,neg_swapped)

sex<-read.table("TSR_sigSNPS_frompermutation2.txt",header = F)
names(sex)<-c("itt",hdr)
sig_neg<-sex[sex$beta<0,]

#do the same for permuted SNPs
neg_effects<-semi_join(ihs_exp,sig_neg,by=c("ps"))
pos_effects<-anti_join(ihs_exp,sig_neg,by=c("ps"))
neg_swapped<-neg_effects[,c(1,2,3,5,4,6,7,8,9,10)]
exp<-rbind(pos_effects,neg_swapped)

test<-obs %>% separate(snp, c("scaf", "pos"), "_")
obs$snp<-test$scaf
nonscaf<-obs[obs$snp != "Scaf5",]

#add classes for plotting
obs$Class<-"Observed"
exp$Class<-"Permuted"
ihs_gw$Class<-"Genome-wide"
nonscaf$Class<-"Non-Scaf5 Obs."
library(reshape2)

#merge
both<-rbind(obs,exp)
all<-rbind(ihs_gw,both)
all1<-rbind(nonscaf,all)
all2<-all1[,c(1,2,3,4,5,6,11)]
long<-melt(all2,id.vars=c("snp","ps","freq","Class"))

#get color palette
library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
#three<-c("#889da8","#b5e9eb","#fcd392",bay[1],bay[3],bay[6])
three<-c(bay[1],bay[3],bay[6],bay[8])

#relevel
long$variable <- factor(long$variable, levels=c("iHS0","iHS1","unstand"))
long$Class <- factor(long$Class, levels=c("Genome-wide","Permuted","Observed","Non-Scaf5 Obs."))

#plot
long %>% filter(variable=="unstand") %>% ggplot(aes(value)) +
  geom_density(alpha=.40, aes( y=..scaled..,color=Class,fill=Class)) +
  theme_classic() +
  scale_color_manual(values = three) +
  scale_fill_manual(values = three) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Density",x="unstandardized iHS") 

long %>% filter(variable=="iHS1") %>% ggplot(aes(value)) +
  geom_density(alpha=.40, aes( y=..scaled..,color=Class,fill=Class)) +
  theme_classic() +
  scale_color_manual(values = three) +
  scale_fill_manual(values = three) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Density",x="standardized iHS") 

ggplot(ihs_obs) +
  geom_point(aes(freq,iHS1))

############################
#get CIs from resampling
############################

x=nrow(ihs_gw)
y=nrow(obs)


perm_resample<-list(); perm_mean1<-list();  perm_mean0<-list();
for (i in 1:1000) {
  interm<-ihs_gw[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_mean1[[i]] <- median(na.omit(interm$iHS1))
  perm_mean0[[i]] <- median(na.omit(interm$unstand))
  
}

gperm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_mean1))
perm_stats$perm_mean2<-unlist(perm_mean0)

names(perm_stats)<-c("mean_iHS","mean_unstand")

quantile(perm_stats$mean_iHS,probs = c(0.05, 0.95))
median(ihs_obs$iHS1)

quantile(perm_stats$mean_unstand,probs = c(0.05, 0.95))
median(ihs_obs$unstand)

