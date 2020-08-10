library(data.table)
library(tidyverse)
filelist <- dir(pattern = "*.out$")

bin <- 5000

for (k in 1:length(filelist)){
        ehh<-fread(filelist[k])
        df2 <- ehh %>%
                mutate(window = V1 %/% bin) %>%
                group_by(window) %>%
                summarise(median_V3 = median(V3,), median_V4 = median(V4),n_snps=n()) %>%
                mutate(start = (window*bin)+1, stop=(window+1)*bin)
        df2$snp<-filelist[k]
        write.table(df2,"EHH_5kwindowed_SLIM_gen10000.txt",append=T,quote=F,row.names=F,col.names=F)
}
