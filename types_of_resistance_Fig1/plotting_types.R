
#initialize map
install.packages("ggmap")
library(ggmap)
water<-get_stamenmap(bbox= c(left = -97, bottom = 35, 
                             right = -76, top = 45), maptype='toner-lite', zoom=6) #Midwest US
water2<-get_stamenmap(bbox= c(left = -83.5, bottom = 41.8, 
                             right = -81, top = 43), maptype='toner-lite', zoom=9) #Southwestern Ontario

library(scatterpie)
library(wesanderson)

p<- ggmap(water) +
  coord_quickmap()
  #geom_point(data= owgps, aes(x = Longitude, y = Latitude, size = 2), alpha=.5, colour="darkred") +
  #geom_path(data=mdat,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #geom_path(data=ca.provinces,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #scale_color_manual(values = c(prop = "#FF0000", suseptible = "#F2AD00", nat="#00A08A")) 

p2<- ggmap(water2) +
  coord_quickmap()
  #geom_point(data= owgps, aes(x = Longitude, y = Latitude, size = 2), alpha=.5, colour="darkred") +
  #geom_path(data=mdat,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #geom_path(data=ca.provinces,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #scale_color_manual(values = c(prop = "#FF0000", suseptible = "#F2AD00", nat="#00A08A")) 

res_phenos<-read.table("pheno_bypop.txt",header=F)
head(res_phenos)
res_phenos_sumbypop<-res_phenos %>% mutate(tsr_summ=case_when(V5==1 & V4 < 2 ~ 1, #just TSR SNP (V5 = 1 status, V4 < 2 refers to copy number < 2)
                                                                      V5==0 & V4 > 2 ~ 2, #just copy num
                                                                      V5==1 & V4 > 2 ~ 3, #both TSR SNP and copy num
                                                                      V5==0 & V4 < 2 & V3 > 1 ~ 4, #NTSR (neither, but resistant (i,e, V3 > 1)
                                                                      V5==0 & V4 < 2 & V3 < 2 ~ 0)) %>% #sus, no mechanism
  group_by(V2,tsr_summ) %>% summarise(n = n())

names(res_phenos_sumbypop)<-c("Pop","State","Count")
library(tidyr)
pop_props_wide <- spread(res_phenos_sumbypop, State, Count)
pop_props_wide[is.na(pop_props_wide)]<-0
pop_props_wide$sum_n<-rowSums(pop_props_wide[,2:6],na.rm = T)
names(pop_props_wide)<-c("Pop","Sus","SNP","Amp","BothTSR","NTSR","samplesize")
pop_props_wide_coord<-inner_join(pop_props_wide,inds,by="Pop")
head(pop_props_wide_coord)
tail(pop_props_wide_coord)

install.packages("scatterpie")
library(scatterpie)
library(PNWColors)
star<-pnw_palette(name="Starfish",n=7,type="continuous")
als574<-c("#59629B","#203096","#6a2fd6","#a779fc","#adb7f7")
als653<-c("#2a663a","#7da888")
ppo<-c("#E69B99","#bd9291","#8f3b39")


names(pop_props_wide_coord)
names(pop_props_wide_coord)<-c("Pop","574_1","574_2","574_3","574_4","574_5","653_6","653_7","PPO_8","PPO_9","PPO_10","sum_n","lat","long","Long2","Lat2")
pop_props_wide_coord[is.na(pop_props_wide_coord)] <- 0

pop_props_wide_coord$Long2<-jitter(pop_props_wide_coord$long,factor = 60)
pop_props_wide_coord$Lat2<-jitter(pop_props_wide_coord$lat,factor = 60)

pop_props_wide_coord$radius <- 6 * abs(rnorm(pop_props_wide_coord$sum_n))

p2 + geom_scatterpie(data=pop_props_wide_coord[1:11,],aes(y = Lat2, x = Long2, group=Pop, r=sqrt(sum_n)/40),
                    cols=c("574_1","574_2","574_3","574_4","574_5","653_6","653_7","PPO_8","PPO_9","PPO_10"),alpha=.9) +
  #coord_equal() +
  #geom_scatterpie_legend(pop_props_wide_coord$radius, x=-80, y=44) +
  geom_scatterpie_legend((pop_props_wide_coord$sum_n[1:11]/100), n=5, 
                         labeller= function(x) round(x*100),
                         x=-81.35, y=42) +
  scale_fill_manual(values=c(als574[1:4],als653,ppo))

square<-function(x) x * 35

p + geom_scatterpie(data=pop_props_wide_coord[12:19,],aes(y = lat, x = long, group=Pop, r=sum_n/35),
                     cols=c("574_1","574_2","574_3","574_4","574_5","653_6","653_7","PPO_8","PPO_9","PPO_10"),alpha=.9) +
  #coord_equal() +
geom_scatterpie_legend((pop_props_wide_coord$sum_n[12:19]/35), n=6, 
                         labeller= function(x) round(x*35),
                         x=-96, y=43.5) +
  scale_fill_manual(values=c(als574[c(1:3,5)],als653,ppo),name="Mutational\nOrigin")


#########
pop_props_wide_coord[,2:6] <- lapply(pop_props_wide_coord[,2:6], as.numeric)
pop_props_wide_coord$Long2<-jitter(pop_props_wide_coord$long,factor = 90)
pop_props_wide_coord$Lat2<-jitter(pop_props_wide_coord$lat,factor = 90)

ont<-pop_props_wide_coord[c(1:11,20:22),]
head(ont)
states<-pop_props_wide_coord[12:19,]
library(viridis)

p2 + geom_scatterpie(data=pop_props_wide_coord,aes(y = Lat2 , x = Long2, group=Pop, r=samplesize/80),
                     cols=c("SNP","Amp","BothTSR","NTSR","Sus"),alpha=.9) +
  scale_fill_viridis(discrete=TRUE) +
  geom_scatterpie_legend((pop_props_wide_coord$samplesize/80), n=5, 
                         labeller= function(x) round(x*100),
                         x=-81.35, y=42)

p + geom_scatterpie(data=states,aes(y = Lat2 , x = Long2, group=Pop, r=samplesize/20),
                     cols=c("SNP","Amp","BothTSR","NTSR","Sus"),alpha=.9) +
  scale_fill_viridis(discrete=TRUE) +
  geom_scatterpie_legend((states$samplesize/20), n=5, 
                         labeller= function(x) round(x*20),
                         x=-96, y=43.5)
