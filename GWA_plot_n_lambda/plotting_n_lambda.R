install.packages("data.table")
install.packages("dplyr")
library(data.table)
library(ggplot2)
library(dplyr)
library(qqman)


sex<-fread("~/forgwas_nodups_tryagain_updated_indsexcluded_TSR.assoc.txt",header=T)
sex$FDR_corr<-p.adjust(sex$p_lrt,method = "fdr")
sex$bon<-p.adjust(sex$p_lrt,method = "bonferroni")


library(tidyverse)
# Theoretical quantile function for -log10(p)
qlog10 <- function(p.values) {
  theoretical <- rank(p.values)/length(p.values)
  return(-log10(theoretical))
}

sex$qlogp<-qlog10(sex$p_lrt)
sex<-sex[order(-qlogp),]

# QQplot with base 10
ggplot(sex, aes(x = qlogp, y = -log10(p_lrt))) + 
  geom_path() +
  geom_abline(intercept = 0, slope = 1)


## calculates lambda
chisq<-qchisq(1-sex$p_lrt,1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq) / sqrt( length(chisq) ) ) )
lambda + (SE.median / qchisq(0.5, 1))
lambda - (SE.median / qchisq(0.5, 1))

#plot
sex %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_point(alpha=.5) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Position in Mb") +
  theme(axis.text.x = element_text(angle = 90))

sex %>% filter(chr==5) %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_point(alpha=.5) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 6.39597,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.260492,lty="dashed",lwd=.25,color="grey") 

sex %>% filter(chr==5) %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_rect(xmin=23.5,xmax=30,ymin=0,ymax=15, alpha=.5, fill="lightgrey",color="lightgrey") +
  geom_point(alpha=.5) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 5.957412,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.241488,lty="dashed",lwd=.25,color="darkgrey") 

sex %>% filter(chr==5) %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_point(alpha=.5) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none") +
  geom_hline(yintercept = 5.957412,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.241488,lty="dashed",lwd=.25,color="grey") 


library(lemon)
grid_arrange_shared_legend(p1,p2, nrow=2, ncol=1, position="right")


