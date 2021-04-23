################################
#Plotting EHH by haplotype and class
################################
library(tidytree)
library(data.table)
library(tidyr)

obs<-fread("EHH_5kwindowed_poppermuted_observed.txt") #output from from EHH windowing script
obs<-obs %>% separate(V7, c("A", "B","chr","ps"))
obs$chr<-gsub('Scaf', '', obs$chr)
obs$chr<-as.integer(as.character(obs$chr))
obs$ps<-as.integer(as.character(obs$ps))
obs$B<-gsub('rep', '', obs$B)

names(obs)<-c("win","ehh1","ehh0","snps","start","end","gen","rep","chr","ps")
se <- function(x) sqrt(var(x)/length(x))

#keep only windows 10Kb on either side of focal allele (i.e. < 11 windows away)
filted<-obs %>% filter(win < 11 & win > -11) %>% group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0))
ehh1<-filted[,c(1,2,4)]
ehh1$class<-"EHH1" #alternate haplotype
ehh0<-filted[,c(1,3,5)]
ehh0$class<-"EHH0" #reference haplotype

#read in full GWA to filter by signficance
sex<-fread("~/nonbinary_twocovs_n155.assoc.txt",header=T)
sex$FDR_corr<-p.adjust(sex$p_lrt,method = "fdr")
sig<-sex[sex$FDR_corr<0.05,]
sig_neg<-sig[sig$beta<0,]

#extract signficant negative/pos effect observed SNPs
neg_effects<-semi_join(obs,sig_neg,by=c("chr","ps"))
pos_effects<-anti_join(obs,sig_neg,by=c("chr","ps"))
neg_swapped<-neg_effects[,c(1,3,2,4,5,6,7,8,9,10)]
obs<-rbind(pos_effects,neg_swapped)

#expected (randomized) pop-permuted SNPs, same deal as above
sex<-read.table("twocovs_FDR05_sigsnps.txt",header = F)
names(sex)<-c("itt",hdr)
sig_neg<-sex[sex$beta<0,]

exp<-fread("EHH_5kwindowed_2cov_poppermuted.txt")
exp<-exp %>% separate(V7, c("A", "B","chr","ps"))
exp$chr<-gsub('Scaf', '', exp$chr)
exp$chr<-as.integer(as.character(exp$chr))
exp$ps<-as.integer(as.character(exp$ps))

neg_effects<-semi_join(exp,sig_neg,by=c("chr","ps"))
pos_effects<-anti_join(exp,sig_neg,by=c("chr","ps"))
neg_swapped<-neg_effects[,c(1,3,2,4,5,6,7,8,9,10)]
exp<-rbind(pos_effects,neg_swapped)
names(exp)<-c("win","ehh1","ehh0","snps","start","end","gen","rep","chr","ps")

#expected (randomized) region-permuted SNPs, same deal as above
sex<-read.table("twocovs_regional_FDR05_sigsnps.txt",header = F)
names(sex)<-c("itt",hdr)
sig_neg<-sex[sex$beta<0,]

exp_reg<-fread("EHH_5kwindowed_2cov_regionpermuted.txt")
exp_reg<-exp_reg %>% separate(V7, c("A", "B","chr","ps"))
exp_reg$chr<-gsub('Scaf', '', exp_reg$chr)
exp_reg$chr<-as.integer(as.character(exp_reg$chr))
exp_reg$ps<-as.integer(as.character(exp_reg$ps))

neg_effects<-semi_join(exp_reg,sig_neg,by=c("chr","ps"))
pos_effects<-anti_join(exp_reg,sig_neg,by=c("chr","ps"))
neg_swapped<-neg_effects[,c(1,3,2,4,5,6,7,8,9,10)]
exp_reg<-rbind(pos_effects,neg_swapped)
names(exp_reg)<-c("win","ehh1","ehh0","snps","start","end","gen","rep","chr","ps")

se <- function(x) sqrt(var(x)/length(x))

filted<-obs %>%  group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0),n_win=n())
exp2<- exp %>% group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0),n_win=n())
exp_reg<- exp_reg %>% group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0),n_win=n())

#keep windows with at least 10 observations
filted<-filted %>% filter(n_win > 10)
exp2<-exp2 %>% filter(n_win > 10)
exp_reg<-exp_reg %>% filter(n_win > 10)

#add classes for plotting
exp_reg$Class<-"Region Permuted"
exp2$Class<-"Pop Permuted"
filted$Class<-"Observed"
#nonscaf$Class<-"Non-Scaf5 Obs."

both<-rbind(exp2,filted)
all<-rbind(exp_reg,both)

ehh1<-all[,c(1,2,4,6,7)] #extract alternate (R) haplotype
ehh2<-all[,c(1,3,5,6,7)] #extract reference (S) haplotype
names(ehh1)<-c("wind","ehh","se","n_win","class")
names(ehh2)<-c("wind","ehh","se","n_win","class")
ehh1$allele<-"R"
ehh2$allele<-"S"

ehh1_2<-rbind(ehh1,ehh2) #long format
ehh1_2$upper<-ehh1_2$ehh+ehh1_2$se
ehh1_2$lower<-ehh1_2$ehh-ehh1_2$se

#long$Class <- factor(ehh1_2$class, levels=c("Genome-wide","Permuted","Observed","Non-Scaf5 Obs."))

library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
three<-c("#d69387","#ebc488","#90d2d4",bay[8],bay[6],bay[3])
ehh1_2$allele <- factor(ehh1_2$allele, levels=c("S","R"))

ehh1_2 %>% 
  ggplot(aes(wind,ehh,color=interaction(class,allele))) +
  geom_point() +
  geom_line(stat="smooth",method = "loess", span = 0.1,
            size = 1.25,
            alpha = 0.8) +
  geom_errorbar(ymax=ehh1_2$upper,ymin=ehh1_2$lower,width=1.2,cex=1,alpha=0.9) +
  theme_classic() +
  #ylim(0,.3) +
  scale_color_manual(values = three) +
  xlab("Distance in Kb from focal site") +
  ylab("EHH") +
  xlim(-10,10)
  
ehh1_2 %>% filter(allele == "R") %>%
  ggplot(aes(upper,color=interaction(class,allele))) +
  geom_density() +
  #geom_line(stat="smooth",method = "loess", span = 0.1,
  #          size = 1.5,
  #          alpha = 0.5) +
  #geom_errorbar(ymax=ehh1_2$upper,ymin=ehh1_2$lower,width=1) +
  theme_classic() +
  #ylim(0,.8) +
  scale_color_manual(values = three) 
  #xlab("Distance in Kb from focal site") +
  #ylab("EHH")
  #xlim(-10,10)

ehh1_2 %>% filter(allele == "R") %>% group_by(class) %>% summarise(mean_upper=mean(upper))
ehh1_2 %>% filter(allele == "R") %>% group_by(class) %>% summarise(mean_upper=mean(lower))


##############
#get resampled stats
###############
head(obs)

obs_test<- obs[obs$win < 11 & obs$win > -11,]
exp_test<- exp[exp$win < 11 & exp$win > -11,]


#%>% filter(win < 11 & win > -11 ) 
obs_avg<-obs_test %>% group_by(chr,ps) %>% summarise(EHH1_mean=median(ehh1),EHH2_mean=median(ehh0),EHH1_SE=se(ehh1),EHH2_SE=se(ehh0))
exp_avg<- exp_test %>% group_by(chr,ps) %>% summarise(EHH1_mean=median(ehh1),EHH2_mean=median(ehh0),EHH1_SE=se(ehh1),EHH2_SE=se(ehh0))

x=nrow(exp_avg)
y=nrow(obs_avg)

perm_resample<-list(); perm_mean1<-list();  perm_mean0<-list();
for (i in 1:1000) {
  interm<-obs_avg[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_mean1[[i]] <- median(interm$EHH1_mean,na.rm = T)
  
}

#perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_mean1))
perm_stats$perm_mean2<-unlist(perm_mean0)

names(perm_stats)<-c("mean_EHH1","mean_EHH2")

quantile(perm_stats$mean_EHH1,probs = c(0.05, 0.95),na.rm = T) #expected 95% CI for EHH1 (res)
median(obs_avg$EHH1_mean) #observed median

quantile(perm_stats$mean_EHH2,probs = c(0.05, 0.95)) #expected 95% CI for EHH1 (sus)
median(obs_avg$EHH2_mean) #observed median 
