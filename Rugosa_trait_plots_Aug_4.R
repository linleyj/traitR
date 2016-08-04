setwd("C:/Users/Clarrien/Documents/traitR/rugosa_traitR")
rug <-read.csv("C:/Users/Clarrien/Documents/traitR/rugosa_traitR/rugosaTraits.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
library(ggplot2)
dev.off()
rug2<-subset(rug, rug$Full=="TRUE")
d201314<-rug2
d1314_alt<-d201314[d201314$Species=="A",]
table(d1314_alt$Site,d1314_alt$Eurosta)
d1314_alt$Site <-factor(d1314_alt$Site, levels = c("SS01", "SS02","SS03", "SS04", "SS05","ME01","ME02","ME03","ME04","ME05","ME06","ME07","ME08","ME09","ME10","NY11","NY12","NY13","NY14","NY15"))

#grand means for each region and trait
alt<-subset(rug2, rug2$Species=="A")
p1.means<-aggregate(stem.pub.length~region2+Eurosta, alt, mean)
p2.means<-aggregate(stem.height~region2+Eurosta, alt, mean)
p3.means<-aggregate(stem.width~region2+Eurosta, alt, mean)
p4.means<-aggregate(leaf.pub.dens~region2+Eurosta, alt, mean)
p5.means<-aggregate(stem.pub.dens~region2+Eurosta, alt, mean)

(p1<-ggplot(d1314_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.pub.length))
+geom_boxplot()+geom_point(aes(colour=Eurosta,drop=F),na.rm=F)+theme_bw()+
scale_x_discrete(name="Site",labels=c
("NB01","","NB02","","NB03","", "NB04","","NB05","","ME01","ME05","","ME06","","ME07","","ME08","","ME09","ME10","NY11","","NY12","NY14",""))
+scale_y_continuous(name="Stem pubescence length")+scale_colour_discrete(name="",breaks=c("a", "una"),
labels=c("Attacked", "Unattacked"))) 
p1+ theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_point(size=3,aes(colour=Eurosta,drop=F))


(p2<-ggplot(d1314_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.height))
+geom_boxplot()+geom_point(aes(colour=Eurosta,drop=F),na.rm=F)+theme_bw()+
  scale_x_discrete(name="Site",labels=c
                   ("SS01","","SS02","","SS03","", "SS04","","SS05","","ME01","ME05","","ME06","","ME07","","ME08","","ME09","ME10","NY11","","NY12","NY14",""))
+scale_y_continuous(name="Ramet height")+scale_colour_discrete(name="",breaks=c("a", "una"),
                                                                         labels=c("Attacked", "Unattacked")))   
p2 + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_point(size=3,aes(colour=Eurosta,drop=F))

(p3<-ggplot(d1314_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.width))
+geom_boxplot()+geom_point(aes(colour=Eurosta,drop=F),na.rm=F)+theme_bw()+
  scale_x_discrete(name="Site",labels=c
                   ("SS01","","SS02","","SS03","", "SS04","","SS05","","ME01","ME05","","ME06","","ME07","","ME08","","ME09","ME10","NY11","","NY12","NY14",""))
+scale_y_continuous(name="Stem width")+scale_colour_discrete(name="",breaks=c("a", "una"),
                                                                         labels=c("Attacked", "Unattacked")))   
p3 + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_point(size=3,aes(colour=Eurosta,drop=F))

(p4<-ggplot(d1314_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=leaf.pub.dens))
+geom_boxplot()+geom_point(aes(colour=Eurosta,drop=F),na.rm=F)+theme_bw()+
  scale_x_discrete(name="Site",labels=c
                   ("SS01","","SS02","","SS03","", "SS04","","SS05","","ME01","ME05","","ME06","","ME07","","ME08","","ME09","ME10","NY11","","NY12","NY14",""))
+scale_y_continuous(name="Leaf pubescence density")+scale_colour_discrete(name="",breaks=c("a", "una"),
                                                                         labels=c("Attacked", "Unattacked")))   
p4 + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_point(size=3,aes(colour=Eurosta,drop=F))

(p5<-ggplot(d1314_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.pub.dens))
+geom_boxplot()+geom_point(aes(colour=Eurosta,drop=F),na.rm=F)+theme_bw()+
  scale_x_discrete(name="Site",labels=c
                   ("SS01","","SS02","","SS03","", "SS04","","SS05","","ME01","ME05","","ME06","","ME07","","ME08","","ME09","ME10","NY11","","NY12","NY14",""))
+scale_y_continuous(name="Stem pubescence density")+scale_colour_discrete(name="",breaks=c("a", "una"),
                                                                          labels=c("Attacked", "Unattacked")))   
p5 + theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_point(size=3,aes(colour=Eurosta,drop=F))


####################################

full<-subset(d1314_alt, Full=="TRUE")
att<-subset(full,Eurosta=="a")
una<-subset(full,Eurosta=="una")

(a<-mean(att$stem.height)) #105.6817
(a<-mean(una$stem.height)) #107.8442

(b<-mean(att$stem.width)) #3.942
(b<-mean(una$stem.width))  #3.680603

(c<-mean(att$stem.pub.length)) #1.232
(c<-mean(una$stem.pub.length)) #1.183103

(d<-mean(att$leaf.pub.dens)) #32.43833
(d<-mean(una$leaf.pub.dens)) #34.22915

#t-test values by trait

#stem pub length
#2014 New Brunswick
t.test(altOV$stem.pub.length~altOV$Eurosta) #t = 0.26715, df = 80.717,  p-value = 0.79,  
Ua=1.097714  Uuna=1.076531 %diff #1.97
#2013 Maine Central
t.test(altOV$stem.pub.length~altOV$Eurosta)#t = 0.047452, df = 23.717, p-value = 0.9626, 
Ua=1.191176 Uuna=1.186644  %diff #0.38
#2013 All Maine
t.test(altOV$stem.pub.length~altOV$Eurosta)#t = 0.20937, df = 19.238, p-value = 0.8364, 
Ua=1.191176  Uuna=1.172244 %diff #1.62
#2013 New York
t.test(altOV$stem.pub.length~altOV$Eurosta) #t = 1.8385, df = 8.5648, p-value = 0.1008, U
a=1.906250  Uuna=1.470109 %diff #29.67
library(metap)
p<-c(0.79,0.9626,0.8364,0.1008)
(NB<-sqrt(137))
(NY<-sqrt(27))
(ME<-sqrt(90))
(NBa<-sqrt(35))
(NYa<-sqrt(8))
(MEa<-sqrt(17))
w<-c(NB,ME,NY)
w<-c(NBa,MEa,NYa)
x<-sumz(p, weights = w, data = NULL, subset = NULL, na.action = na.fail)
print(x)

install.packages("BSDA")
library(BSDA)
data=c(55, 58, 61, 61, 62, 62, 62, 63, 63, 64, 66, 68, 68, 69, 69, 69, 70, 71, 72, 72)
sign.test(data, conf.level=0.90)

#ramet height
#2014 NB
t.test(altOV$stem.height~altOV$Eurosta) #t = -0.90028, df = 66.99, p-value = 0.3712,  
((97.29429/101.34694)-1)*100 %diff #-3.998789
#2013 Maine all
t.test(altOV$stem.height~altOV$Eurosta) #t = 0.87216, df = 24.909, p-value = 0.3915,  
((113.4353/109.8260)-1)*100  %diff #3.28638
#2013 NY
t.test(altOV$stem.height~altOV$Eurosta) #t = 1.4215, df = 12.691, p-value = 0.1793, 
((125.9000/110.7435)-1)*100  %diff #13.68613

#stem width
#2014 NB
t.test(altOV$stem.width~altOV$Eurosta) #t = -0.56868, df = 70.679, p-value = 0.5714, 
((3.931429/4.075510)-1)*100  %diff #-3.535288
#2013 Maine all
t.test(altOV$stem.width~altOV$Eurosta) #t = 0.66345, df = 23.465, p-value = 0.5135, 
((3.630000/3.515512)-1)*100 %diff #3.256652
#2013 NY
t.test(altOV$stem.width~altOV$Eurosta) #t = 2.2587, df = 23.016, p-value = 0.03368, 
((4.65125/3.75087)-1)*100  %diff #24.00456

#leaf pub dens
#2014 NB
t.test(altOV$leaf.pub.dens~altOV$Eurosta) #t = -0.77629, df = 69.243, p-value = 0.4402, 
((7.237143/8.930612)-1)*100 %diff#-18.96252
#2013 Maine
t.test(altOV$leaf.pub.dens~altOV$Eurosta) #t = 0.65864, df = 18.653, p-value = 0.5182,
((42.58824/33.00787)-1)*100 %diff#29.0245
#2013 NY
t.test(altOV$leaf.pub.dens~altOV$Eurosta) #t = 1.1398, df = 21.412, p-value = 0.267, 
((121.12500/94.86957)-1)*100  %diff#27.67529

#stem pub dens
#2014 NB
t.test(altOV$stem.pub.dens~altOV$Eurosta) #t = 0.0073022, df = 49.52, p-value = 0.9942
((194.4829/194.2878)-1)*100 #0.100418
#2013 Maine
t.test(altOV$stem.pub.dens~altOV$Eurosta) #t = 0.4554, df = 18.865, p-value = 0.654
((160.0840/151.5298)-1)*100 #5.645226
#2013 NY
t.test(altOV$stem.pub.dens~altOV$Eurosta)#t = -1.1907, df = 13.656, p-value = 0.2541
((121.7857/150.9938)-1)*100 #-19.34391

rug <-read.csv("C:/Users/Clarrien/Documents/traitR/rugosa_traitR/rugosaTraits.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
str(rug)
rug2<-subset(rug, rug$Full=="TRUE")
dim(rug2)


#2013 New York
NY<-rug2[rug2$region=="NY",]
dim(NY)
altOV<-NY[NY$Species=="A",]
dim(altOV)
t.test(altOV$stem.pub.length~altOV$Eurosta) #t = 1.8385, df = 8.5648, p-value = 0.1008, 
((1.906250/1.470109)-1)*100 #29.66726
t.test(altOV$leaf.pub.dens~altOV$Eurosta) #t = 1.1398, df = 21.412, p-value = 0.267, 
((121.12500/94.86957)-1)*100  #27.67529
t.test(altOV$stem.height~altOV$Eurosta) #t = 1.4215, df = 12.691, p-value = 0.1793, 
((125.9000/110.7435)-1)*100  #13.68613
t.test(altOV$stem.width~altOV$Eurosta) #t = 2.2587, df = 23.016, p-value = 0.03368, 
((4.65125/3.75087)-1)*100  #24.00456
t.test(altOV$stem.pub.dens~altOV$Eurosta)#t = -1.1907, df = 13.656, p-value = 0.2541
((121.7857/150.9938)-1)*100 #-19.34391
#2013 Maine-east - no attacked alt so no data
MEe<-rug2[rug2$region=="ME-east",]
dim(MEe)
altOV<-MEe[MEe$Species=="A",]
dim(altOV)
str(altOV)
t.test(altOV$stem.pub.length~altOV$Eurosta) 
#2013 Maine-central
MEc<-rug2[rug2$region=="ME-cent",]
dim(MEc)
altOV<-MEc[MEc$Species=="A",]
dim(altOV)
str(altOV)
t.test(altOV$stem.pub.length~altOV$Eurosta) #t = 0.047452, df = 23.717, p-value = 0.9626, 1.191176  1.186644 
#2013 Maine-west - no attacked altissima 
MEw<-rug2[rug2$region=="ME-west",]
dim(MEw)
altOV<-MEw[MEw$Species=="A",]
dim(altOV)
str(altOV)
t.test(altOV$stem.pub.length~altOV$Eurosta)
#2014 SS
NB<-rug2[rug2$region=="SS",]
dim(NB)
altOV<-NB[NB$Species=="A",]
dim(altOV)
str(altOV)
t.test(altOV$stem.pub.length~altOV$Eurosta) #t = 0.26715, df = 80.717, p-value = 0.79, 
((1.097714/1.076531)-1)*100 #1.967709
t.test(altOV$leaf.pub.dens~altOV$Eurosta) #t = -0.77629, df = 69.243, p-value = 0.4402, 
((7.237143/8.930612)-1)*100 #-18.96252
t.test(altOV$stem.height~altOV$Eurosta) #t = -0.90028, df = 66.99, p-value = 0.3712,  
((97.29429/101.34694)-1)*100 #-3.998789
t.test(altOV$stem.width~altOV$Eurosta) #t = -0.56868, df = 70.679, p-value = 0.5714, 
((3.931429/4.075510)-1)*100  #-3.535288
t.test(altOV$stem.pub.dens~altOV$Eurosta) #t = 0.0073022, df = 49.52, p-value = 0.9942
((194.4829/194.2878)-1)*100 #0.100418

#2013 All Maine
ME<-rug2[rug2$region2=="ME",]
dim(ME)
altOV<-ME[ME$Species=="A",]
dim(altOV)
t.test(altOV$stem.pub.length~altOV$Eurosta)#t = 0.20937, df = 19.238, p-value = 0.8364, 
((1.191176/1.172244)-1)*100 #1.615022
t.test(altOV$leaf.pub.dens~altOV$Eurosta) #t = 0.65864, df = 18.653, p-value = 0.5182,
((42.58824/33.00787)-1)*100 #29.0245
t.test(altOV$stem.height~altOV$Eurosta) #t = 0.87216, df = 24.909, p-value = 0.3915,  
((113.4353/109.8260)-1)*100  #3.28638
t.test(altOV$stem.width~altOV$Eurosta) #t = 0.66345, df = 23.465, p-value = 0.5135, 
((3.630000/3.515512)-1)*100 #3.256652
t.test(altOV$stem.pub.dens~altOV$Eurosta) #t = 0.4554, df = 18.865, p-value = 0.654
((160.0840/151.5298)-1)*100 #5.645226


#correlations for 2013+2014 data
library(psych)
alt<-subset(rug2,Species=="A")
str(alt)
Alt.traits<-subset(alt[c(8,9,10,11,16)])
Alt.traits.corr.coefficients<-corr.test(Alt.traits)$r
Alt.traits.corr.pvalues<-corr.test(Alt.traits)$p
write.csv(Alt.traits.corr.coefficients,"201314Rugosa.Alt.corr.coefficients.csv",quote=F,row.names=F)
write.csv(Alt.traits.corr.pvalues,"201314Rugosa.Alt.corr.pvalues.csv",quote=F,row.names=F)

#regional correlations
rug <-read.csv("C:/Users/Clarrien/Documents/traitR/rugosa_traitR/rugosaTraits.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
str(rug)
rug2<-subset(rug, rug$Full=="TRUE")
dim(rug2)

#all maine
ME<-rug2[rug2$region2=="ME",]
dim(ME)
alt<-ME[ME$Species=="A",]
dim(alt)
Alt.traits<-subset(alt[c(8,9,10,11,16)])
Alt.traits.corr.coefficients<-corr.test(Alt.traits)$r
Alt.traits.corr.pvalues<-corr.test(Alt.traits)$p
write.csv(Alt.traits.corr.coefficients,"Maine.Alt.corr.coefficients.csv",quote=F,row.names=F)
write.csv(Alt.traits.corr.pvalues,"Maine.Alt.corr.pvalues.csv",quote=F,row.names=F)

#New Brunswick
NB<-rug2[rug2$region2=="SS",]
dim(NB)
alt<-NB[NB$Species=="A",]
dim(alt)
Alt.traits<-subset(alt[c(8,9,10,11,16)])
Alt.traits.corr.coefficients<-corr.test(Alt.traits)$r
Alt.traits.corr.pvalues<-corr.test(Alt.traits)$p
write.csv(Alt.traits.corr.coefficients,"NB.Alt.corr.coefficients.csv",quote=F,row.names=F)
write.csv(Alt.traits.corr.pvalues,"NB.Alt.corr.pvalues.csv",quote=F,row.names=F)

#NEw York
NY<-rug2[rug2$region2=="NY",]
dim(NY)
alt<-NY[NY$Species=="A",]
dim(alt)
Alt.traits<-subset(alt[c(8,9,10,11,16)])
Alt.traits.corr.coefficients<-corr.test(Alt.traits)$r
Alt.traits.corr.pvalues<-corr.test(Alt.traits)$p
write.csv(Alt.traits.corr.coefficients,"NY.Alt.corr.coefficients.csv",quote=F,row.names=F)
write.csv(Alt.traits.corr.pvalues,"NY.Alt.corr.pvalues.csv",quote=F,row.names=F)

#try weighted tests
install.packages("weights")
library(weights)
rug <-read.csv("C:/Users/Clarrien/Documents/traitR/rugosa_traitR/rugosaTraits.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
rug2<-subset(rug, rug$Full=="TRUE")
alt<-subset(rug2, rug2$Species=="A")
attach(alt)
spl.means<-aggregate(stem.pub.length~Site+Eurosta, alt, mean)
spl.var<-aggregate(stem.pub.length~Site+Eurosta, alt, var)
spl.sd<-aggregate(stem.pub.length~Site+Eurosta, alt, sd)
spl.mean.var.sd<-cbind(spl.means,spl.var$stem.pub.length,spl.sd$stem.pub.length)
spl.mean.var.sd$var.weight<-(1/spl.var$stem.pub.length)
spl.mean.var.sd$sqrt.sd.weight<-(1/(sqrt(spl.sd$stem.pub.length)))
head(spl.mean.var.sd)
write.csv(spl.mean.var.sd,"stem_pub_length_weights.csv",quote=F,row.names=F)

spl.df<-read.csv("C:/Users/Clarrien/Documents/traitR/rugosa_traitR/stem_pub_length_weights.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
spl.df
spl.df2<-spl.df[c(1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,19,21,22,23,24,25,26),]
wtd.t.test(spl.df$mean.diff, weight=spl.df$sqrt.sd.weight)
wtd.t.test(spl.df$mean.diff, weight=spl.df$n)

tryaov<-aov(stem.pub.length~Eurosta,data=spl.df2,weights=var.weight)
summary(tryaov)
#try z test, so need t-tests for every site? goddamit
site<-subset(alt,Site=="NY11") #t = 1.1878, df = 6.161, p-value = 0.2787, 1.812500          1.455357 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="NY14") #t = 1.0464, df = 2.634, p-value = 0.3818, 2.187500          1.729167 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS01")  #t = 0.15909, df = 12.979, p-value = 0.876, 1.0025000         0.9881818 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS02") #t = 1.0521, df = 15.891, p-value = 0.3085, 1.118889          1.024545 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS03") #t = -0.16292, df = 1.2224, p-value = 0.8933, 1.030             1.078 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS04")  #t = 0.79555, df = 11.179, p-value = 0.4429, 1.170000          1.052727 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS05")  #t = -0.61712, df = 15.716, p-value = 0.546, 1.083846          1.240000 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME05")  #t = 0.051479, df = 13.159, p-value = 0.9597, 1.13750           1.13125 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME06") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME07") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME08") #t = 0.49252, df = 5.3288, p-value = 0.642, 1.3375            1.2325 
t.test(site$stem.pub.length~site$Eurosta)
p<-c(0.2787,0.3818,0.876,0.3085,0.8933,0.4429,0.546,0.9597,0.642)
site<-subset(alt,Site=="NY11") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="NY14") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS01") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS02") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS03") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS04") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="SS05") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME05") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME06") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME07") 
t.test(site$stem.pub.length~site$Eurosta)
site<-subset(alt,Site=="ME08") 
t.test(site$stem.pub.length~site$Eurosta)



install.packages("metap")
library(metap)
p<-c(0.2787,0.3818,0.876,0.3085,0.8933,0.4429,0.546,0.9597,0.642)
w<-c(3.16227766, 2.236067977, 2.449489743,1.414213562, 2, 3,1.414213562,2.645751311, 3.605551275)
x<-sumz(p, weights = w, data = NULL, subset = NULL, na.action = na.fail)
print(x)


d<-mean(3.75,5.125,1.1875,5.5,1.125)
d
stem.pub.dens.means<-aggregate(stem.pub.dens~Site+Eurosta, rug2, mean)
leaf.pub.dens.means<-aggregate(leaf.pub.dens~Site+Eurosta, rug2, mean)
stem.height.means<-aggregate(stem.height~Site+Eurosta, rug2, mean)
stem.width.means<-aggregate(stem.width~Site+Eurosta, rug2, mean)

?

una<-subset(rug2,Eurosta=="una")
a<-subset(rug2,Eurosta=="a")



str(rug)
rug$weights<-(rug$stem.pub.length)
2<-subset(rug, rug$Full=="TRUE")
dim(rug2)

#try sign tests
library(coin)
attach(rug2)
sign_test(stem.pub.length~Eurosta,data=rug2,subset=Site=="SS01",alternative="greater")
sign_test(stem.pub.length~Eurosta,data=rug2,subset=Site=="SS02",alternative="greater")
sign_test(stem.pub.length~Eurosta,data=rug2,subset=Site=="SS02")
sign_test(stem.pub.length~Eurosta,data=rug2,subset=Site=="SS03")
sign_test(stem.pub.length~Eurosta,data=rug2,subset=Site=="SS04",alternative="greater")
sign_test(stem.pub.length~Eurosta,data=rug2,subset=Site=="SS05",alternative="greater")

#fisher exact
Input =(
  "Trait	NB	ME	NY
Galled	3	2	2
  ungalled	2	2	0
  ")

Matriz = as.matrix(read.table(textConnection(Input),  header=TRUE,    row.names=1))
Matriz
fisher.test(Matriz, alternative="two.sided")

#try sign tests
spl.allsites<-binom.test(3, 3) #success is higher for galled  
print(spl.allsites)..

rug <-read.csv("C:/Users/Clarrien/Documents/traitR/rugosa_traitR/rugosaTraits.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
str(rug)
library(ggplot2)
dev.off()
d2013<-rug[rug$Year=="2013",]
table(d2013$Site,d2013$Eurosta)
str(d2013)
d13_alt<-d2013[d2013$Species=="A",]

p1<-ggplot(d13_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.pub.length*leaf.pub.dens))+geom_boxplot()+geom_point(aes(colour=Eurosta, drop=F),na.rm=F)+scale_x_discrete("Site",drop=F)+theme_bw()
p1
p2<-ggplot(d13_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=leaf.pub.dens))+geom_boxplot()+geom_point(aes(colour=Eurosta, drop=F),na.rm=F)+scale_x_discrete("Site",drop=F)+theme_bw()

(p3<-ggplot(d13_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.height))+geom_boxplot()+geom_point(aes(colour=Eurosta, drop=F),na.rm=F)+scale_x_discrete("Site",drop=F)+theme_bw()+theme(panel.background = element_blank()) )

(p4<-ggplot(d13_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.width))+geom_boxplot()+geom_point(aes(colour=Eurosta, drop=F),na.rm=F)+scale_x_discrete("Site",drop=F)+theme_bw())
p1
?ggsave
pdf("AltissimaTraitPlots2013.pdf")
print(p1)
p2
p3
p4
dev.off()


rug <-read.csv("C:/Users/Clarrien/Documents/traitR/rugosa_traitR/rugosaTraits.csv",  strip.white = TRUE,    na.strings = c("NA", ""))
rug$stem.pub.length.scalecent<-(rug$stem.pub.length-mean(rug$stem.pub.length,na.rm = T))/sd(rug$stem.pub.length,na.rm = T)
rug$stem.pub.dens.scalecent<-(rug$stem.pub.dens-mean(rug$stem.pub.dens,na.rm = T))/sd(rug$stem.pub.dens,na.rm = T)

rug<-subset(rug, rug$Full=="TRUE")
d201314<-rug
table(d201314$Site,d201314$Eurosta)

str(d201314)
d201314$Site<-as.factor(d201314$Site)
levels(d201314$Site)

d201314$Site <-factor(d201314$Site, levels = c("SS01", "SS02","SS03", "SS04", "SS05","ME01","ME02","ME03","ME04","ME05","ME06","ME07","ME08","ME09","ME10","NY11","NY12","NY13","NY14","NY15"))
d1314_alt<-d201314[d201314$Species=="A",]
head(d1314_alt) 

rug$stem.pub.length.scalecent<-(rug$stem.pub.length-mean(rug$stem.pub.length,na.rm = T))/sd(rug$stem.pub.length,na.rm = T)
rug$stem.pub.dens.scalecent<-(rug$stem.pub.dens-mean(rug$stem.pub.dens,na.rm = T))/sd(rug$stem.pub.dens,na.rm = T)


(p<-ggplot(d1314_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.pub.length.scalecent*stem.pub.dens.scalecent))+geom_boxplot()+geom_point(aes(colour=Eurosta, drop=F),na.rm=F)+scale_x_discrete("Site",drop=F)+theme_bw())

#Linley final code
(p1<-ggplot(d1314_alt,aes(x=interaction(Eurosta,as.factor(Site)),y=stem.pub.length))+geom_boxplot()+geom_point(aes(colour=Eurosta,drop=F),na.rm=F)+theme_bw()+
  scale_x_discrete(name="Site",labels=c("SS01",
                                        "","SS02","","SS03","", "SS04",
                                        "","SS05","","ME05","","ME06","","ME07","","ME08","","ME09","","ME10","","NY11","","NY12","","","NY13","","NY14","","NY15"))+scale_y_continuous(name="Pubescence")+scale_colour_discrete(name="",
                                                                                                                                                                                                                                 breaks=c("a", "una"),
                                                                                                                                                                                                                                 
