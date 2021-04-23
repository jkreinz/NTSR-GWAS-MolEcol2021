library(ggplot2)
library(tidyverse)
library(PNWColors)
library(data.table)

#read in observed, expected (permuted), and genomewide sets
ihs_obs<-fread("allobserved_IHS.txt")
ihs_exp<-fread("2cov_regionalperm_IHSpermutedhits.txt")
ihs_gw<-fread("allscafs.ihs")
ihs_exp_pop<-fread("2cov_popperm_IHSpermutedhits.txt")

names(ihs_obs)<-c("snp","ps","freq","iHS1","iHS0","unstand","left","right","ancestral_left","ancestral_right")
names(ihs_exp)<-c("snp","ps","freq","iHS1","iHS0","unstand","left","right","ancestral_left","ancestral_right")
names(ihs_gw)<-c("snp","ps","freq","iHS1","iHS0","unstand","left","right","ancestral_left","ancestral_right")
names(ihs_exp_pop)<-c("snp","ps","freq","iHS1","iHS0","unstand","left","right","ancestral_left","ancestral_right")

#ihs_obs<-ihs_obs %>% separate(snp, c("chr", "ps"), "_")
#ihs_obs$chr<-as.integer(gsub("Scaf", "", ihs_obs$chr))
#ihs_obs$ps<-as.integer(as.character(ihs_obs$ps))
#sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
#ihs_obs<-left_join(ihs_obs,sex,by=c("chr","ps"))


#swap iHH1 and iHH0 values such that all 1 alleles have positive effects on resistance
sex<-fread("~/nonbinary_twocovs_n155.assoc.txt",header=T)
hdr<-names(sex)
sex$FDR_corr<-p.adjust(sex$p_lrt,method = "fdr")
sig<-sex[sex$FDR_corr<0.05,]
sig_neg<-sig[sig$beta<0,]

neg_effects<-semi_join(ihs_obs,sig_neg,by=c("ps"))
pos_effects<-anti_join(ihs_obs,sig_neg,by=c("ps"))
neg_swapped<-neg_effects[,c(1,2,3,5,4,6,7,8,9,10)]
obs<-rbind(pos_effects,neg_swapped)

sex<-read.table("twocovs_regional_FDR05_sigsnps.txt",header = F)
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

### and for pop permuted set

sex<-read.table("twocovs_FDR05_sigsnps.txt",header = F)
names(sex)<-c("itt",hdr)
sig_neg<-sex[sex$beta<0,]

#do the same for permuted SNPs
neg_effects<-semi_join(ihs_exp_pop,sig_neg,by=c("ps"))
pos_effects<-anti_join(ihs_exp_pop,sig_neg,by=c("ps"))
neg_swapped<-neg_effects[,c(1,2,3,5,4,6,7,8,9,10)]
exp_pop<-rbind(pos_effects,neg_swapped)

#exp_pop$Class<-"Non-Scaf5 Obs."

#######################
test<-obs %>% separate(snp, c("scaf", "pos"), "_")
obs$snp<-test$scaf
nonscaf<-obs[obs$snp != "Scaf5",]

#add classes for plotting
obs$Class<-"Observed"
exp$Class<-"Region Permuted"
ihs_gw$Class<-"Genome-wide"
nonscaf$Class<-"Non-Scaf5 Obs."
exp_pop$Class<-"Pop Permuted"

library(reshape2)

#merge
both<-rbind(obs,exp)
all<-rbind(ihs_gw,both)
all<-rbind(all,exp_pop)
all2<-all[,c(1,2,3,4,5,6,11)]
long<-melt(all2,id.vars=c("snp","ps","freq","Class"))

#STANDARDIZE
stand<-read.table("standardized_ihs_popperm.txt",skip = 1)
names(ihs_obs)
names(stand)<-c("freq","n","mean_ihs","var_ihs")

obs$freq<-ceiling(obs$freq*100)/100
obs_wstand<-inner_join(obs,stand,by = "freq")
obs_wstand$stand<-(obs_wstand$unstand-obs_wstand$mean_ihs)/obs_wstand$var_ihs

exp$freq<-ceiling(exp$freq*100)/100
exp_wstand<-inner_join(exp,stand,by = "freq")
exp_wstand$stand<-(exp_wstand$unstand-exp_wstand$mean_ihs)/exp_wstand$var_ihs

exp_pop$freq<-ceiling(exp_pop$freq*100)/100
exp_pop_wstand<-inner_join(exp_pop,stand,by = "freq")
exp_pop_wstand$stand<-(exp_pop_wstand$unstand-exp_pop_wstand$mean_ihs)/exp_pop_wstand$var_ihs


ihs_gw$freq<-ceiling(ihs_gw$freq*100)/100
ihsgw_wstand<-inner_join(ihs_gw,stand,by = "freq")
ihsgw_wstand$stand<-(ihsgw_wstand$unstand-ihsgw_wstand$mean_ihs)/ihsgw_wstand$var_ihs

test<-exp_wstand %>% separate(snp, c("scaf", "pos"), "_")
exp_wstand$snp<-test$scaf

test<-ihsgw_wstand  %>% separate(snp, c("scaf", "pos"), "_")
ihsgw_wstand$snp<-test$scaf

#add classes for plotting
obs_wstand$Class<-"Observed"
exp_wstand$Class<-"Region Permuted"
exp_pop_wstand$Class<-"Pop Permuted"
ihsgw_wstand$Class<-"Genome-wide"
#nonscaf_wstand$Class<-"Non-Scaf5 Obs."


library(reshape2)
#merge
both<-rbind(obs_wstand,exp_wstand)
all<-rbind(ihsgw_wstand,both)
all<-rbind(all,exp_pop_wstand)

all2<-all[,c(1,2,3,6,11,15)]
long<-melt(all2,id.vars=c("snp","ps","freq","Class"))

IHHs<-all[,c(1,2,3,4,5,11)]
long2<-melt(IHHs,id.vars=c("snp","ps","freq","Class"))


#get color palette
library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
three<-c("#889da8","#b5e9eb","#fcd392","",bay[1],bay[3],bay[6])
three<-c(bay[1],bay[3],bay[6],bay[8])

#relevel
#long$variable <- factor(long$variable, levels=c("iHS0","iHS1","unstand"))
long$Class <- factor(long$Class, levels=c("Genome-wide","Region Permuted","Pop Permuted","Observed"))

#plot
long %>% filter(variable=="unstand") %>% ggplot(aes(value)) +
  geom_density(alpha=.40, aes( y=..scaled..,color=Class,fill=Class)) +
  theme_classic() +
  scale_color_manual(values = three) +
  scale_fill_manual(values = three) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Density",x="unstandardized iHS") 

long %>% filter(variable=="stand") %>% ggplot(aes(value)) +
  geom_density(alpha=.40, aes( y=..scaled..,color=Class,fill=Class)) +
  theme_classic() +
  scale_color_manual(values = three) +
  scale_fill_manual(values = three) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Density",x="standardized iHS") +
  xlim(-5,+5)


long2$variable <- factor(long2$variable, levels=c("iHS1","iHS0"))
long$Class <- factor(long2$Class, levels=c("Genome-wide","Region Permuted","Pop Permuted","Observed"))
three2<-c("#889da8",bay[1],"#b5e9eb",bay[3],"#fcd392",bay[6],"#ebafa4",bay[8])
three2<-c(bay[8],bay[6],bay[3],"#ebafa4","#fcd392","#b8cccc")

long2 %>% filter(Class=="Observed" | Class=="Pop Permuted" | Class=="Region Permuted") %>% ggplot(aes(value)) +
  geom_density(alpha=.4, aes( y=..scaled..,color=interaction(Class,variable),fill=interaction(Class,variable))) +
  theme_classic() +
  scale_color_manual(values = three2) +
  scale_fill_manual(values = three2) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Density",x="iHH") +
  xlim(0,.5)


############################
#get CIs from resampling
############################

(x=nrow(exp_pop_wstand))
(y=nrow(obs_wstand))


perm_resample<-list(); perm_mean1<-list();  perm_mean0<-list();
for (i in 1:1000) {
  interm<-exp_pop_wstand[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_mean1[[i]] <- median(interm$stand,na.rm = T)
  perm_mean0[[i]] <- median(interm$unstand,na.rm = T)
  
}

#gperm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_mean1))
perm_stats$perm_mean2<-unlist(perm_mean0)

names(perm_stats)<-c("mean_stand","mean_unstand")

quantile(perm_stats$mean_stand,probs = c(0.05, 0.95))
median(obs_wstand$stand)
mean(obs_wstand$stand)

quantile(perm_stats$mean_unstand,probs = c(0.05, 0.95))
median(obs_wstand$unstand)


perm_resample<-list(); perm_mean1<-list();  perm_mean0<-list();
for (i in 1:1000) {
  interm<-exp_pop_wstand[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_mean1[[i]] <- median(interm$iHS1,na.rm = T)
  perm_mean0[[i]] <- median(interm$iHS0,na.rm = T)
  
}

gperm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_mean1))
perm_stats$perm_mean2<-unlist(perm_mean0)

names(perm_stats)<-c("iHS1","iHS0")

quantile(perm_stats$iHS1,probs = c(0.05, 0.95))
median(obs_wstand$iHS1,na.rm = T)

quantile(perm_stats$iHS0,probs = c(0.05, 0.95))
median(obs_wstand$unstand)
