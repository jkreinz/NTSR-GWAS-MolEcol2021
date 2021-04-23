#for no covariate model
for i in {1..250}; do awk -v var="$i" '($14 < 0.000000402) {print var "\t" $0}' forgwas_nodups_tryagain_updated_indsexcluded_vanilla_${i}.assoc.txt; done >> vanilla_sigSNPS_frompermutation.txt &
#value is equivalent to pvalue cutoff at q < 0.05

#for +TSR model
for i in {1..250}; do awk -v var="$i" '($14 < 0.000001102) {print var "\t" $0}' forgwas_nodups_tryagain_updated_indsexcluded_TSR_${i}.assoc.txt; done >> TSR_sigSNPS_frompermutation.txt &
#value is equivalent to pvalue cutoff at q < 0.05

#for +TSR+EPSPS model
for i in {1..250}; do awk -v var="$i" '($14 < 0.00000251188) {print var "\t" $0}' forgwas_nodups_tryagain_updated_indsexcluded_twocovs_${i}.assoc.txt; done >> twocovs_sigSNPS_frompermutation.txt &
#value is equivalent to pvalue cutoff at q < 0.05

########
#read into R and calculate number of average # significant hits per GWAS (estimated FDR)
R
TSR_sigs<-read.table("vanilla_sigSNPS_frompermutation.txt") #all significant hits across all permuted GWASs
hdr<-read.table("forgwas_nodups_tryagain_updated_indsexcluded_TSR_1.assoc.txt",header=T,nrows=1) #header
names(TSR_sigs)<-c("itt",names(hdr))

library(dplyr)
counts_byitt<-table(TSR_sigs$itt) #table of number significant hits @ q < 0.05 across GWAS iterations

ittcounts<-as.vector(counts_byitt)
ittcounts<-c(ittcounts,rep(0,(250-length(ittcounts)))) #add 0's for however many GWAS itts had no sig hits

average_sighits_perGWA<-mean(ittcounts) #get mean FDR
median_sighits_perGWA<-median(ittcounts) # get median FDR

#######################

