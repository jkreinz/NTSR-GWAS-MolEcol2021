for i in {1..250}; do awk -v var="$i" '($14 < 0.000000402) {print var "\t" $0}' forgwas_nodups_tryagain_updated_indsexcluded_vanilla_${i}.assoc.txt; done >> vanilla_sigSNPS_frompermutation.txt &
#value is equivalent to pvalue cutoff at q < 0.05

for i in {1..250}; do awk -v var="$i" '($14 < 0.000001102) {print var "\t" $0}' forgwas_nodups_tryagain_updated_indsexcluded_TSR_${i}.assoc.txt; done >> TSR_sigSNPS_frompermutation.txt &
#value is equivalent to pvalue cutoff at q < 0.05


########
#calculate number of significant hits per GWAS
R
TSR_sigs<-read.table("vanilla_sigSNPS_frompermutation.txt")
hdr<-read.table("forgwas_nodups_tryagain_updated_indsexcluded_TSR_1.assoc.txt",header=T,nrows=1)
names(TSR_sigs)<-c("itt",names(hdr))

library(dplyr)
counts_byitt<-table(TSR_sigs$itt)
TSR_meannumsig<-sum(counts_byitt)/250

ittcounts<-as.vector(counts_byitt)
ittcounts<-c(ittcounts,rep(0,(250-length(ittcounts))))

average_sighits_perGWA<-mean(ittcounts)
median_sighits_perGWA<-median(ittcounts)

#######################

