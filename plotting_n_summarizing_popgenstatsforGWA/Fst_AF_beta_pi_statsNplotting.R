
library(data.table)
library(dplyr)

#read in between R and S Fst values
fst<-fread("forgwas_nodups_tryagain_updated_indsexcluded_RvS.weir.fst",header=T)
fst[fst$WEIR_AND_COCKERHAM_FST<0,]$WEIR_AND_COCKERHAM_FST<-0
hist(fst$WEIR_AND_COCKERHAM_FST)
mean(na.omit(fst$WEIR_AND_COCKERHAM_FST))

#read in gff, retain genes with sig hits
hitsbygene<-read.table("ntsr_sig_tsr_n155_sorted.gff",fill = T)
hitsbygene<-hitsbygene[,1:13]
hitsbygene<-hitsbygene[complete.cases(hitsbygene),]
hitsbygene <- hitsbygene[hitsbygene$V10 == ".",]

hitsbygene$V1<-gsub("[a-zA-Z_]", "", hitsbygene$V1)
names(hitsbygene)[c(1,3)]<-c("CHROM","POS")
hitsbygene$CHROM<-as.integer(hitsbygene$CHROM)
hitsbygene$CHROM<-as.integer(hitsbygene$CHROM)
hitsbygene$POS<-as.integer(hitsbygene$POS)

#merge gene names of sig hits and fst values
TSRsig<-semi_join(fst,hitsbygene,by=c("CHROM","POS"))
TSRsig_not5<-TSRsig[TSRsig$CHROM!=5,]
TSRsig_not<-anti_join(fst,hitsbygene,by=c("CHROM","POS"))

#classes for plotting
TSRsig_not$Class<-"Insig."
TSRsig_not5$Class<-"Non-scaf5 Sig."
TSRsig$Class<-"Sig."

#permuted fst values
perm_fst<-fread("TSRallsig_fakephenos.fst",header=T)
names(perm_fst)<-names(fst)
perm_fst$Class<-"Permuted"

#merge
both<-rbind(TSRsig_not,perm_fst)
all<-rbind(TSRsig,both)
more<-rbind(TSRsig_not5,all)

library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[1],bay[3],bay[6],bay[8])

library(ggplot2)

more$Class <- factor(more$Class, levels=c("Insig.","Permuted","Sig.","Non-scaf5 Sig."))
ggplot(data=more, aes(WEIR_AND_COCKERHAM_FST)) +
  geom_density(aes(y=..scaled.., fill=Class,color=Class), alpha=.4) +
  theme_classic() +
  xlab("Weir & Cockerham Fst") +
  ylab("Frequency") +
  scale_fill_manual(values=three) +
  scale_color_manual(values=three)
  

#get resampled 95% stats

x <- nrow(perm_fst)
y <- 196

perm_fst[perm_fst$WEIR_AND_COCKERHAM_FST<0,]$WEIR_AND_COCKERHAM_FST<-0

perm_resample<-list(); perm_mean<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-perm_fst[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_mean[[i]] <- median(na.omit(interm$WEIR_AND_COCKERHAM_FST))
  perm_variance[[i]] <- var(na.omit(interm$WEIR_AND_COCKERHAM_FST))
  
}

perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_mean))
names(perm_stats)<-"mean_fst"
perm_stats$variance<-unlist(perm_variance)

quantile(perm_stats$mean_fst,probs = c(0.05, 0.95))
median(TSRsig$WEIR_AND_COCKERHAM_FST)


#check association between AF and fst
perm_fst<-fread("TSRallsig_fakephenos.fst",header=F)
perm_AF<-read.table("TSR_sigSNPS_frompermutation2.txt",header = F)
names(perm_AF)[c(2,4,8)]<-c("X1","X2","af")
names(perm_fst)[c(1,2,3)]<-c("X1","X2","fst")

permed<-cbind(perm_fst,perm_AF[,c(2,4,8)])

sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
names(sex)[c(1,3)]<-c("X1","X2")
sex<-sex[,c(1,3,7)]
TSRsig_AF<-semi_join(sex,hitsbygene,by=c("CHROM","POS"))
observed<-cbind(TSRsig,TSRsig_AF)

observed$class<-"Sig"
permed$class<-"Permuted"
both<-rbind(observed,permed,use.names=F)
names(both)<-c("Chrom","Pos","Fst","X1","X2","AF","Class")

bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[3],bay[6])


ggplot(both, aes(AF,Fst, color=Class)) +
  geom_point(alpha=.2) +
  geom_smooth() +
  scale_color_manual(values=three) +
  theme_classic()

########################
#AF
########################


library(data.table)
library(dplyr)
library(ggplot2)

sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
#sex$i<-500
#sex$FDR_corr<-p.adjust(sex$p_lrt,method = "fdr")
names(sex)[c(1,3)]<-c("CHROM","POS")
sex<-sex[,c(1,3,7)]

permuted_sigSNPS<-read.table("TSR_sigSNPS_frompermutation2.txt",header = F)
names(permuted_sigSNPS)[c(2,4,8)]<-c("CHROM","POS","af")
permuted_sigSNPS<-permuted_sigSNPS[,c(2,4,8)]

TSRsig<-semi_join(sex,hitsbygene,by=c("CHROM","POS"))
TSRsig_not5<-TSRsig[TSRsig$CHROM!=5,]
TSRsig_not<-anti_join(sex,hitsbygene,by=c("CHROM","POS"))


TSRsig_not$Class<-"Insig."
TSRsig_not5$Class<-"Non-scaf5 Sig."
TSRsig$Class<-"Sig."
permuted_sigSNPS$Class<-"Permuted"

both<-rbind(TSRsig_not,permuted_sigSNPS)
all<-rbind(TSRsig,both)
more<-rbind(TSRsig_not5,all)

library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[1],bay[3],bay[6],bay[8])

#library(ggplot2)

more$Class <- factor(more$Class, levels=c("Insig.","Permuted","Sig.","Non-scaf5 Sig."))

ggplot(data=more, aes(af)) +
  geom_density(alpha=.5, aes(y=..scaled..,color=Class,fill=Class)) +
  theme_classic() +
  scale_color_manual(values = three) +
  scale_fill_manual(values = three) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y="Density",x="Allele frequency")

#get stats
x=nrow(permuted_sigSNPS)
y=196

perm_resample<-list(); perm_mean<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-permuted_sigSNPS[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_mean[[i]] <- median(na.omit(interm$af))
  perm_variance[[i]] <- var(na.omit(interm$af))
  
}

perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_mean))
names(perm_stats)<-"mean_AF"
perm_stats$variance<-unlist(perm_variance)

quantile(perm_stats$mean_AF,probs = c(0.05, 0.95))
median(TSRsig$af)

#################
#PI
#################


fst<-fread("forgwas_nodups_tryagain_updated_indsexcluded_RES.sites.pi",header=T)
hist(btw$PI)
mean(na.omit(fst$PI))

TSRsig<-semi_join(fst,hitsbygene,by=c("CHROM","POS"))
TSRsig_not<-anti_join(fst,hitsbygene,by=c("CHROM","POS"))
TSRsig_not$Class<-"Insig."
TSRsig$Class<-"Sig."
TSRsig_not5<-TSRsig[TSRsig$CHROM!=5,]
TSRsig_not5$Class<-"Non-scaf5 Sig."

perm_fst<-fread("TSR_permutedsigsites.sites.pi",header=T)
perm_fst$Class<-"Permuted"

both<-rbind(TSRsig_not,TSRsig)
all<-rbind(perm_fst,both)
more<-rbind(TSRsig_not5,all)

bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[1],bay[3],bay[6],bay[8])

more$Class <- factor(more$Class, levels=c("Insig.","Permuted","Sig.","Non-scaf5 Sig."))

ggplot(data=more, aes(PI)) +
  geom_density(aes(y=..scaled.., fill=Class,color=Class), alpha=.5) +
  theme_classic() +
  xlab("SNP-wise Nucleotide Diversity") +
  ylab("Frequency") +
  scale_fill_manual(values=three) +
  scale_color_manual(values=three)
#geom_vline(xintercept = median(TSRsig$WEIR_AND_COCKERHAM_FST), lty="dashed",cex=1) +
#geom_vline(xintercept = quantile(na.omit(fst$WEIR_AND_COCKERHAM_FST),probs = (0.999)), color="grey", lty="dashed")

#get stats

x <- nrow(perm_fst)
y <- 196

perm_resample<-list(); perm_median<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-perm_fst[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_median[[i]] <- median(na.omit(interm$PI))
  perm_variance[[i]] <- var(na.omit(interm$PI))
  
}

perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_median))
names(perm_stats)<-"median_pi"
perm_stats$variance<-unlist(perm_variance)

quantile(perm_stats$median_pi,probs = c(0.05, 0.95))
median(TSRsig$PI)
################
#pi between
################

btw<-fread("forgwas_nodups_tryagain_updated_indsexcluded_allinds.sites.pi",header=T)
res_fst<-fread("forgwas_nodups_tryagain_updated_indsexcluded_res.sites.pi",header=T)
sus_fst<-fread("forgwas_nodups_tryagain_updated_indsexcluded_sus.sites.pi",header=T)

within_avg<-(res_fst$PI+sus_fst$PI)/2
btw$pi_boverw<-btw$PI/within_avg

perm_bwn<-fread("sus_res_all_allreps.sites.pi",header=F)
names(perm_bwn)<-c("CHROM","POS","Sus","chr2","pos2","Res","chr3","pos3","All")
perm_bwn$avg_within<-(perm_bwn$Sus+perm_bwn$Res)/2

perm_bwn$pi_boverw<-perm_bwn$All/perm_bwn$avg_within
perm_bwn$Class<-"Permuted"


TSRsig<-semi_join(btw,hitsbygene,by=c("CHROM","POS"))
TSRsig_not<-anti_join(btw,hitsbygene,by=c("CHROM","POS"))
TSRsig_not$Class<-"Insig."
TSRsig$Class<-"Sig."
TSRsig_not5<-TSRsig[TSRsig$CHROM!=5,]
TSRsig_not5$Class<-"Non-scaf5 Sig."

perm_bwn<-perm_bwn[,c(1,2,6,11,12)]

both<-rbind(TSRsig_not,TSRsig)
all<-rbind(perm_bwn,both,use.names=F)
more<-rbind(TSRsig_not5,all,use.names=F)

bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[1],bay[3],bay[6],bay[8])

more$Class <- factor(more$Class, levels=c("Insig.","Permuted","Sig.","Non-scaf5 Sig."))

ggplot(data=more, aes(pi_boverw)) +
  geom_density(aes(y=..scaled.., fill=Class,color=Class), alpha=.5) +
  theme_classic() +
  xlab("SNP-wise Pi Between/Pi Within") +
  ylab("Frequency") +
  scale_fill_manual(values=three) +
  scale_color_manual(values=three) 



 #get stats

x <- nrow(perm_bwn)
y <- 196

perm_resample<-list(); perm_median<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-perm_bwn[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  perm_median[[i]] <- median(na.omit(interm$pi_boverw))
  perm_variance[[i]] <- var(na.omit(interm$pi_boverw))
  
}

perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_median))
names(perm_stats)<-"median_pib"
perm_stats$variance<-unlist(perm_variance)

quantile(perm_stats$median_pib,probs = c(0.05, 0.95))

######################
#BETA

sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
#sex$FDR_corr<-p.adjust(sex$p_lrt,method = "fdr")
#sig<-sex[sex$FDR_corr<0.05,]

names(sex)[c(1,3)]<-c("CHROM","POS")

permuted_sigSNPS<-read.table("TSR_sigSNPS_frompermutation2.txt",header = F)
names(permuted_sigSNPS)[c(2,4,8)]<-c("CHROM","POS","af")

TSRsig<-semi_join(sex,hitsbygene,by=c("CHROM","POS"))
TSRsig_not5<-TSRsig[TSRsig$CHROM!=5,]
TSRsig_not<-anti_join(sex,hitsbygene,by=c("CHROM","POS"))


#TSRsig_not$Class<-"Insig."
TSRsig_not5$Class<-"Non-scaf5 Obs."
TSRsig$Class<-"Observed"
TSRsig<-TSRsig[,c(1,3,7,8,9,16)]
TSRsig_not5<-TSRsig_not5[,c(1,3,7,8,9,16)]

permuted_sigSNPS$Class<-"Permuted"
permuted_sigSNPS<-permuted_sigSNPS[,c(2,4,8,9,10,17)]

names(permuted_sigSNPS)<-names(TSRsig)
both<-rbind(permuted_sigSNPS,TSRsig)
all<-rbind(both,TSRsig_not5)

all$upper<-all$beta+all$se
all$lower<-all$beta-all$se

bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[3],bay[6],bay[8])

all$Class <- factor(all$Class, levels=c("Permuted","Observed","Non-scaf5 Obs."))

ggplot(data=all, aes(af,beta,fill=Class,color=Class)) +
  geom_errorbar(ymax=all$upper,ymin=all$lower,alpha=.2) +
  geom_point(alpha=.4,cex=4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=three) +
  scale_color_manual(values=three) +
  labs(y="Effect size (Beta)",x="Allele frequency")


x <- nrow(permuted_sigSNPS)
y <- 196

perm_resample<-list(); perm_median<-list(); perm_variance<-list()
for (i in 1:1000) {
  interm<-permuted_sigSNPS[sample(x,y,replace = T),]
  perm_resample[[i]] <- interm
  neg<-sum(interm$beta<0)
  pos<-sum(interm$beta>=0)
  perm_median[[i]] <-neg/(pos+neg)

}

perm_resample_df<-do.call(rbind, perm_resample)
perm_stats<-as.data.frame(do.call(rbind, perm_median))
names(perm_stats)<-"propneg_beta"
#perm_stats$variance<-unlist(perm_variance)

quantile(perm_stats$propneg_beta,probs = c(0.05, 0.95))
neg<-sum(TSRsig$beta<0)
pos<-sum(TSRsig$beta>=0)
neg/(pos+neg)

x <- nrow(permuted_sigSNPS)
y <- 196

####
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

quantile(perm_stats$median_beta_abs,probs = c(0.05, 0.95))


median(abs(TSRsig_not5$beta))


####################
#get number of significant hits by gene
####################

hitsbygene<-read.table("ntsr_sig_tsr_n155_sorted.gff",fill = T)
hitsbygene<-hitsbygene[,1:13]
hitsbygene<-hitsbygene[complete.cases(hitsbygene),]

gwas_TSR<-hitsbygene %>% group_by(V13) %>% tally
#gwas_TSR$V1 <- factor(gwas_TSR$V1, levels = c("Scaffold_1","Scaffold_2","Scaffold_3","Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"))

ggplot(data=gwas_TSR) +
  geom_histogram(aes(n))

ggplot(data=gwas_TSR, aes(n)) +
  geom_histogram(fill="grey", color="black",binwidth = 1) +
  theme_classic() +
  xlab("Significant hits per gene") +
  ylab("Frequency")


library(ggplot2)

hitsbygene_nocov<-read.table("NTSR_vanilla_n155.gff",fill = T, na.strings = "")
hitsbygene_nocov<-hitsbygene_nocov[,1:13]
hitsbygene_nocov<-hitsbygene_nocov[complete.cases(hitsbygene_nocov),]
hitsbygene_nocov <- hitsbygene_nocov[hitsbygene_nocov$V10 == ".",]

test<-hitsbygene_nocov %>% group_by(V13) %>% tally
hist(test$n, main="", xlab="Significant hits per gene", xlim=c(1,7))


gwas_nocov<-hitsbygene_nocov %>% group_by(V1) %>% tally
gwas_nocov$V1 <- factor(gwas_nocov$V1, levels = c("Scaffold_1","Scaffold_2","Scaffold_3","Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"))


ggplot(data=gwas_nocov, aes(V1,n)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))  

ggplot(data=gwas_TSR, aes(V1,n)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))  

gwas_TSR$model<-"TSR"
gwas_TSR$prop<-gwas_TSR$n/(sum(gwas_TSR$n))

gwas_nocov$model<-"no_covariate"
gwas_nocov$prop<-gwas_nocov$n/(sum(gwas_nocov$n))

both<-rbind(gwas_TSR,gwas_nocov)

hitsbygene_CN<-read.table("gwas_wCN_outliers_sorted.gff",fill = T, na.strings = "")
hitsbygene_CN<-hitsbygene_CN[,1:13]
hitsbygene_CN<-hitsbygene_CN[complete.cases(hitsbygene_CN),]
hitsbygene_CN <- hitsbygene_CN[hitsbygene_CN$V10 == ".",]

gwas_CN<-hitsbygene_CN %>% group_by(V1) %>% tally
gwas_CN$model<-"copy_number"
gwas_CN$V1 <- factor(gwas_CN$V1, levels = c("Scaffold_1","Scaffold_2","Scaffold_3","Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"))

three<-rbind(both,gwas_CN)


hitsbygene_both<-read.table("NTSR_bothcov_n155.ggg",fill = T, na.strings = "")
hitsbygene_both<-hitsbygene_both[,1:13]
hitsbygene_both<-hitsbygene_both[complete.cases(hitsbygene_both),]
hitsbygene_both <- hitsbygene_both[hitsbygene_both$V10 == ".",]

gwas_both<-hitsbygene_both %>% group_by(V1) %>% tally
gwas_both$model<-"both"
gwas_both$V1 <- factor(gwas_both$V1, levels = c("Scaffold_1","Scaffold_2","Scaffold_3","Scaffold_4","Scaffold_5","Scaffold_6","Scaffold_7","Scaffold_8","Scaffold_9","Scaffold_10","Scaffold_11","Scaffold_12","Scaffold_13","Scaffold_14","Scaffold_15","Scaffold_16"))


four<-rbind(three,gwas_both)

install.packages("PNWColors")
library(PNWColors)
bay<-pnw_palette("Bay",8,type="continuous")
three<-c(bay[4],bay[2])

ggplot(data=both, aes(V1,prop, color=model)) +
  geom_point(cex=3,alpha=.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))  +
  labs(x="",y="Proportion of Sig. SNPs") +
  scale_color_manual(values=three)



ggplot(data=both, aes(V1,n, color=model)) +
  geom_point(cex=2,alpha=.8) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="",y="Frequency")

