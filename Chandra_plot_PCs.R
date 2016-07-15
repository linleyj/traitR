installed.packages("devEMF")
library(emf)
update.packages("ggplot2")
library(grid)
library(ggplot2)
library(dplyr)
library(reshape2)
d1<-read.csv("~/Downloads/Gnor.June.PCoA.multi3_GNOR.csv")
d2<-read.csv("~/Downloads/Gnor.phen.numeric.multi3.noForce_GNOR.csv")
d3<-read.csv("~/Downloads/June.PCoA.multi3_redone.csv")
d4<-read.csv("~/Downloads/Late_PCA_3axis_herbs.csv")
d5<-read.csv("~/Downloads/May.PCA.multi.all.csv")
head(d1)
d1$Time<-rep("Early",dim(d1)[1])
d1$variable<-rep("GNOR",dim(d1)[1])
names(d1)[2]<-"value"
d2$Time<-rep("Late",dim(d2)[1])
d2$variable<-rep("GNOR",dim(d2)[1])
names(d2)[2]<-"value"
s3<-select(d3,-starts_with("Species."))
head(s3)
m3<-melt(s3,id.vars = c("Species","Comp.1","Comp.2","Comp.3"),measure.vars=5:11)
m3$Time<-rep("Early",dim(m3)[1])
m3$Species.Attacked<-interaction(m3$Species,m3$value)
levels(m3$Species.Attacked)

s4<-select(d4,-starts_with("Species."))
head(s4)
m4<-melt(s4,id.vars = c("Species","Comp.1","Comp.2","Comp.3"),measure.vars=5:11)
m4$Time<-rep("Late",dim(m4)[1])
m4$Species.Attacked<-interaction(m4$Species,m4$value)
levels(m4$Species.Attacked)

head(d5)
s5<-select(d5,-starts_with("Species."))
head(s5)
m5<-melt(s5,id.vars = c("Species","Comp.1","Comp.2"),measure.vars=5:11)
m5$Comp.3<-rep(1,dim(m5)[1])
m5$Species.Attacked<-interaction(m5$Species,m5$value)
m5$Time<-rep("Initial",dim(m5)[1])

bigdat<-rbind(d1,d2,m3,m4,m5)
first3<-bigdat %>%
filter(variable=="EUR"|variable=="GNOR"|variable=="RHOP")
first3$Time<-factor(first3$Time,levels=c("Initial","Early","Late"))
first3$variable<-factor(first3$variable,levels=c("EUR","GNOR","RHOP"))

second4<-bigdat %>%
filter(variable=="Asphondylia.solidaginis"|variable=="Asteromyia"|variable=="Exema"|variable=="Nemorimyza")
second4$Time<-factor(second4$Time,levels=c("Initial","Early","Late"))
second4$variable<-factor(second4$variable)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plotfun<-function(data)
{
  ggplot(data,aes(x=Comp.1,y=Comp.2))+
  #add points scaled by Comp.3
  geom_point(aes(col=Species.Attacked),shape=16,size=0.5)+
    scale_fill_manual(cbPalette)+
  scale_size(range=c(0.01,0.5))+
  #facets split over time and bug
  facet_wrap(~Time+variable, scales ="free")+
  #black and white theme
  theme_classic()+
  #with modifications
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),strip.background = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks = element_blank(),strip.text=element_blank(),axis.title.y=element_text(margin=margin(0,20,0,0)))+
 #add the contours
   geom_density_2d(data=subset(data,Species.Attacked=="A.a"),alpha=0.5,colour=cbPalette[1],show.legend = F)+
    geom_density_2d(data=subset(data,Species.Attacked=="A.una"),alpha=0.5,colour="#E69F00",show.legend = F)+
    geom_density_2d(data=subset(data,Species.Attacked=="G.a"),alpha=0.5,colour=cbPalette[3],show.legend = F)+
   geom_density_2d(data=subset(data,Species.Attacked=="G.una"),alpha=0.5,colour=cbPalette[4],show.legend = F)+
  #change the names of the legend
  scale_color_discrete(breaks=c("A.a",  "A.una","G.a","G.una"),labels=c("Altissima attacked",  "Altissima unattacked","Rugosa attacked", "Rugosa unattacked"))
  }

#set up a plotting device with a constant size.
install.packages("ReporteRs")
library(ReporteRs)


p1<-plotfun(first3)
print(p1)
x.axis.limits = ggplot_build(p1)$panel$ranges[[1]][["x.range"]]
y.axis.limits = ggplot_build(p1)$panel$ranges[[1]][["y.range"]]
pushViewport(dataViewport(yscale = c(0,1), xscale = c(0,1), clip = "off"))
grid.text("Eurosta", x = 0.18, y = 0.98, just = "center", gp = gpar(fontsize=8,col = "black",fontface="italic"), default.units = "native")
grid.text("Gnorimoschema", x = 0.4, y = 0.98, just = "center", gp = gpar(fontsize=8,col = "black",fontface="italic"), default.units = "native")
grid.text("Rhopalomyia", x = 0.6, y = 0.98, just = "center", gp = gpar(fontsize=8,col = "black",fontface="italic"), default.units = "native")
grid.text("Early", x = 0.07, y = 0.8, just = "center", gp = gpar(fontsize=8,col = "black"), default.units = "native",rot=90)
grid.text("Mid", x = 0.07, y = 0.5, just = "center", gp = gpar(fontsize=8,col = "black"), default.units = "native",rot=90)
grid.text("Late", x = 0.07, y = 0.2, just = "center", gp = gpar(fontsize=8,col = "black"), default.units = "native",rot=90)

#this doesn't seem to keep the labels
#but exports an editable graph
mydoc = pptx(  )
mydoc = addSlide( mydoc, slide.layout = "Title and Content" )
myplot<-last_plot()
mydoc = addPlot( mydoc, function( ) print( myplot ), vector.graphic=TRUE)
writeDoc( mydoc, file = "test plot.pptx" )

p1<-plotfun(second4)
print(p1)
x.axis.limits = ggplot_build(p1)$panel$ranges[[1]][["x.range"]]
y.axis.limits = ggplot_build(p1)$panel$ranges[[1]][["y.range"]]
pushViewport(dataViewport(yscale = c(0,1), xscale = c(0,1), clip = "off"))
grid.text("Asphondylia", x = 0.18, y = 0.98, just = "center", gp = gpar(fontsize=8,col = "black",fontface="italic"), default.units = "native")
grid.text("Asteromyia", x = 0.33, y = 0.98, just = "center", gp = gpar(fontsize=8,col = "black",fontface="italic"), default.units = "native")
grid.text("Exema", x = 0.46, y = 0.98, just = "center", gp = gpar(fontsize=8,col = "black",fontface="italic"), default.units = "native")
grid.text("Nemorimyza", x = 0.6, y = 0.98, just = "center", gp = gpar(fontsize=8,col = "black",fontface="italic"), default.units = "native")
grid.text("Early", x = 0.07, y = 0.8, just = "center", gp = gpar(fontsize=8,col = "black"), default.units = "native",rot=90)
grid.text("Mid", x = 0.07, y = 0.5, just = "center", gp = gpar(fontsize=8,col = "black"), default.units = "native",rot=90)
grid.text("Late", x = 0.07, y = 0.2, just = "center", gp = gpar(fontsize=8,col = "black"), default.units = "native",rot=90)
dev.off()
