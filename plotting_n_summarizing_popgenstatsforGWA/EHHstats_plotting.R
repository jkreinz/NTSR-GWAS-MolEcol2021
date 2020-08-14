################################
################################

library(data.table)
obs<-fread("EHH_5kwindowed_allTSRobserved.txt") #read in 5000bp window EHHs

obs<-obs %>% separate(V7, c("A", "B","chr","ps")) #get scafs output SNP name
obs$chr<-gsub('Scaf', '', obs$chr)
obs$chr<-as.integer(as.character(obs$chr))
obs$ps<-as.integer(as.character(obs$ps))

obs$B<-gsub('rep', '', obs$B)
names(obs)<-c("win","ehh1","ehh0","snps","start","end","gen","rep","chr","ps") #name columns
se <- function(x) sqrt(var(x)/length(x)) #define SE

#get windows 100kb around focal snps
filted<-obs %>% filter(win < 11 & win > -11) %>% group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0))

#change to long format
ehh1<-filted[,c(1,2,4)]
ehh1$class<-"EHH1"
ehh0<-filted[,c(1,3,5)]
ehh0$class<-"EHH0"

#read in GWAS output
sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
sex$FDR_corr<-p.adjust(sex$p_lrt,method = "fdr")
sig<-sex[sex$FDR_corr<0.05,]
sig_neg<-sig[sig$beta<0,]

#swap 0,1 alleles for significant hit if negative effect on resistance 
neg_effects<-semi_join(obs,sig_neg,by=c("chr","ps"))
pos_effects<-anti_join(obs,sig_neg,by=c("chr","ps"))
neg_swapped<-neg_effects[,c(1,3,2,4,5,6,7,8,9,10)]
obs<-rbind(pos_effects,neg_swapped)

###do the sanme for permuted
sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
hdr<-names(sex)
sex<-read.table("TSR_sigSNPS_frompermutation2.txt",header = F)
names(sex)<-c("itt",hdr)
sig_neg<-sex[sex$beta<0,]

#EHH values from empirical null sig hits
exp<-fread("EHH_5kwindowed_allTSRpermuted2.txt")
exp<-exp %>% separate(V7, c("A", "B","chr","ps"))
exp$chr<-gsub('Scaf', '', exp$chr)
exp$chr<-as.integer(as.character(exp$chr))
exp$ps<-as.integer(as.character(exp$ps))

#swap
neg_effects<-semi_join(exp,sig_neg,by=c("chr","ps"))
pos_effects<-anti_join(exp,sig_neg,by=c("chr","ps"))
neg_swapped<-neg_effects[,c(1,3,2,4,5,6,7,8,9,10)]
exp<-rbind(pos_effects,neg_swapped)
names(exp)<-c("win","ehh1","ehh0","snps","start","end","gen","rep","chr","ps")

se <- function(x) sqrt(var(x)/length(x))

#get windows
filted<-obs %>%  group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0),n_win=n())
exp2<- exp %>% group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0),n_win=n())
nonscaf<-obs %>%  filter(gen != "Scaf5") %>% group_by(win) %>% summarise(EHH1_mean=mean(ehh1),EHH0_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH0_SE=se(ehh0),n_win=n())

filted<-filted %>% filter(n_win > 39)
exp2<-exp2 %>% filter(n_win > 39)
nonscaf<-nonscaf %>% filter(n_win > 39)

exp2$Class<-"Permuted"
filted$Class<-"Observed"
nonscaf$Class<-"Non-Scaf5 Obs."

#join datasets for plotting
both<-rbind(exp2,filted)
all<-rbind(nonscaf,both)

ehh1<-all[,c(1,2,4,6)]
ehh2<-all[,c(1,3,5,6)]
names(ehh1)<-c("wind","ehh","se","n_win","class")
names(ehh2)<-c("wind","ehh","se","n_win","class")
ehh1$allele<-"R"
ehh2$allele<-"S"

ehh1_2<-rbind(ehh1,ehh2)
ehh1_2$upper<-ehh1_2$ehh+ehh1_2$se
ehh1_2$lower<-ehh1_2$ehh-ehh1_2$se

#long$Class <- factor(ehh1_2$class, levels=c("Genome-wide","Permuted","Observed","Non-Scaf5 Obs."))
#color scheme
library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[8],bay[6],bay[3],"#deafa6","#fcd392","#b5e9eb")

#plot
ehh1_2 %>% 
  ggplot(aes(wind,ehh,color=interaction(class,allele))) +
  geom_point() +
  #geom_line(stat="smooth",method = "loess", span = 0.1,
  #          size = 1.5,
  #          alpha = 0.5) +
  geom_errorbar(ymax=ehh1_2$upper,ymin=ehh1_2$lower,width=1) +
  theme_classic() +
  ylim(0,.45) +
  scale_color_manual(values = three) +
  xlab("Distance in Kb from focal site") +
  ylab("EHH") 
  

##############
#get resampled stats
###############
head(obs)

#obs_test<- obs[obs$win < -1 | obs$win > 1,] #check 1kb around
#exp_test<- exp[exp$win < -1 | exp$win > 1,] #check 1kb around

#%>% filter(win < 11 & win > -11 ) 
obs_avg<-obs_test %>% filter(win < 11 & win > -11 ) %>% group_by(chr,ps) %>% summarise(EHH1_mean=mean(ehh1),EHH2_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH2_SE=se(ehh0))
exp_avg<- exp_test %>% filter(win < 11 & win > -11 ) %>% group_by(chr,ps) %>% summarise(EHH1_mean=mean(ehh1),EHH2_mean=mean(ehh0),EHH1_SE=se(ehh1),EHH2_SE=se(ehh0))
x=nrow(exp_avg)
y=nrow(obs_avg)

#resample w replacement to n(obs.GWAS) 1000 times
perm_resample<-list(); perm_mean1<-list();  perm_mean0<-list();
for (i in 1:1000) {
  interm<-exp_avg[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_mean1[[i]] <- median(na.omit(interm$EHH1_mean))
  perm_mean0[[i]] <- median(na.omit(interm$EHH2_mean))
  
}

perm_resample_df<-do.call(rbind, perm_resample) 
perm_stats<-as.data.frame(do.call(rbind, perm_mean1))
perm_stats$perm_mean2<-unlist(perm_mean0)

names(perm_stats)<-c("mean_EHH1","mean_EHH2")

quantile(perm_stats$mean_EHH1,probs = c(0.05, 0.95)) #median EHH1 95% CIs
median(obs_avg$EHH1_mean) #median EHH1 obs

quantile(perm_stats$mean_EHH2,probs = c(0.05, 0.95)) #median EHH2 95% CIs
median(obs_avg$EHH2_mean)  #median EHH2 obs

