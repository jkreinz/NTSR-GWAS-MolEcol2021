install.packages("data.table")
install.packages("dplyr")
library(data.table)
library(ggplot2)
library(dplyr)
library(qqman)

library(tidyverse)
# Theoretical quantile function for -log10(p)
qlog10 <- function(p.values) {
  theoretical <- rank(p.values)/length(p.values)
  return(-log10(theoretical))
}

#read in GEMMA output
res<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_twocovs.assoc.txt",header=T) #two covariate GEMMA GWA results
res$FDR_corr<-p.adjust(res$p_lrt,method = "fdr")
res$bon<-p.adjust(res$p_lrt,method = "bonferroni")

res$qlogp<-qlog10(res$p_lrt)
res<-res[order(-qlogp),]

#read in GEMMA output
res_vanilla<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_vanilla.assoc.txt",header=T) #no covariate GEMMA GWA results
res_vanilla$FDR_corr<-p.adjust(res$p_lrt,method = "fdr")
res_vanilla$bon<-p.adjust(res$p_lrt,method = "bonferroni")

res_vanilla$qlogp<-qlog10(res_vanilla$p_lrt)
res_vanilla<-res_vanilla[order(-qlogp),]

# QQplot with base 10
ggplot(res, aes(x = qlogp, y = -log10(p_lrt))) + 
  geom_path() +
  geom_abline(intercept = 0, slope = 1)


## calculates lambda
chisq<-qchisq(1-res$p_lrt,1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq) / sqrt( length(chisq) ) ) )
lambda + (SE.median / qchisq(0.5, 1))
lambda - (SE.median / qchisq(0.5, 1))

#plot
p1<- res_vanilla %>%
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_point(alpha=.5) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 6.39597,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.260492,lty="dashed",lwd=.25,color="grey")  #plot GW manhattan 

res_vanilla %>% filter(chr==5) %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_rect(xmin=23.5,xmax=30,ymin=0,ymax=15, alpha=.5, fill="lightgrey",color="lightgrey") +
  geom_point(alpha=.5) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 6.39597,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.260492,lty="dashed",lwd=.25,color="grey")  #just chromosome 5, outlining EPSPS amplification

p2<- res %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_point(alpha=.5) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none") +
  geom_hline(yintercept = 5.6,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.241488,lty="dashed",lwd=.25,color="grey") 

res %>% filter(chr==1) %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_rect(xmin=23.5,xmax=30,ymin=0,ymax=15, alpha=.5, fill="lightgrey",color="lightgrey") +
  geom_point(alpha=.5) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 6.39597,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.260492,lty="dashed",lwd=.25,color="grey")  #just chromosome 1, where high prop of sig SNPs map

#to plot no covariate and TSR manhattans together
library(lemon)
grid_arrange_shared_legend(p1,p2, nrow=2, ncol=1, position="right")


