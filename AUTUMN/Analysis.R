##############################GENERAL INTRO##############################


### Set your working directory
setwd(dir="C:\\Users\\darth\\Desktop\\REASEARCH DATA\\AUTUMN")
### Remove all the objects
rm(list=ls(all=TRUE)) 
### Save a backup
#save.image("Backup.RData")
#load("Backup.RData")


### Install and load required packages


#install.packages("iNEXT",dependencies=T)
library(iNEXT)
#install.packages("ggplot2",dependencies=T)
library(ggplot2)
#install.packages("plotrix",dependencies=T)
library(plotrix)
#install.packages("car",dependencies=T)
library(car)
#install.packages("AER",dependencies=T)
library(AER)
#install.packages("hillR",dependencies=T)
library(hillR)
#install.packages("performance",dependencies=T)
library(performance)
#install.packages("sf",dependencies=T)
library(sf)
#install.packages("spdep",dependencies=T)
library(spdep)


### Load dataset


abb <-read.table("Max-Abundance.txt")
env <-read.table("ENV.txt")
env$Category <-as.factor(env$Category)
clc <-read.table("CLC.txt")
clc$Category <-factor(clc$Category)


##############################EXTRAPOLATION - COVERAGE##############################


### The first step is to make the abundance data collected in the plots comparable
### We also get an idea of sample completeness thanks to coverage

DataInfo(abb,datatype="abundance")
set.seed(123)
abbcove <-estimateD(abb,q=c(0,1,2),datatype="abundance",base="coverage",nboot=1000)

### We can also visualize the differences for every plot or for a whole group
#abbrho <-abb[,c(-3:-24)]
#set.seed(123)
#rho <-iNEXT(abbrho,q=c(0,1,2),endpoint=168.95601,datatype="abundance")
#ggiNEXT(rho, type=3, se=TRUE, facet.var="Order.q", color.var="Assemblage")+xlim(0,0.8686516) 
#abbwet <-abb[,c(-2,-4,-6,-8,-10,-12,-14,-16,-18,-20,-22,-24)]
#set.seed(123)
#wet <-iNEXT(abbwet,q=0,endpoint=168.95601,datatype="abundance")
#ggiNEXT(wet, type=3, se=TRUE, color.var="Assemblage")+xlim(0,0.8686516)+theme(legend.position = "none")


##############################MORAN TEST##############################


shp <- st_read("C:\\Users\\darth\\Desktop\\MATTEO\\RICERCA\\Treatment Wetlands\\QGIS Autunno\\Plot 300x300\\Plot 300x300.shp")
shp #not the same order of abbcove
shp_ord <- shp[c(11,12,10,9,8,7,6,5,16,15,14,13,23,24,22,21,17,18,19,20,4,3,2,1), ]
shp_ord
shp_centroids <- st_centroid(shp_ord)
coords <- st_coordinates(shp_centroids)
coords
nb <- dnearneigh(coords, d1 = 0, d2 = 5000)
lw <- nb2listw(nb, style = "W")
moran.test(abbcove[abbcove$Order.q==0,]$qD, lw, alternative = "two.sided", rank=T)
moran.test(abbcove[abbcove$Order.q==1,]$qD, lw, alternative = "two.sided", rank=F)
moran.test(abbcove[abbcove$Order.q==2,]$qD, lw, alternative = "two.sided", rank=F)

shp_treat <-shp[c(11,10,8,6,16,14,23,22,17,19,4,2), ]
shp_centrotreat <- st_centroid(shp_treat)
coordstreat <- st_coordinates(shp_centrotreat)
nbt <- dnearneigh(coordstreat, d1 = 0, d2 = 5000)
lwt <- nb2listw(nbt, style = "W", zero.policy=T)
moran.test(covewetS, lwt, alternative = "two.sided", rank=T)
moran.test(covewetH, lwt, alternative = "two.sided", rank=F)
moran.test(covewetD, lwt, alternative = "two.sided", rank=F)


##############################PAIRED TEST - WETLANDS VS CONTROL##############################


### Coverage-based Richness

covewetS <-abbcove[abbcove$Order.q==0,]$qD[c(1,3,5,7,9,11,13,15,17,19,21,23)]
covecontrS <-abbcove[abbcove$Order.q==0,]$qD[c(2,4,6,8,10,12,14,16,18,20,22,24)]
wilcox.test(covewetS,covecontrS,paired=T)
boxplot(covewetS,covecontrS)
coveS <-abbcove[abbcove$Order.q==0,]$qD
cove05S <-abbcove[abbcove$Order.q==0,]$qD.LCL
cove95S <-abbcove[abbcove$Order.q==0,]$qD.UCL
par(mar = c(2, 5, 2, 2) + 0.1)
plotCI(x=c(1:24),y=coveS,li=cove05S,ui=cove95S,lwd=2,cex.lab=2,cex.axis=1.5,pch=21,scol=c("#D81B60","#1E88E5"),xaxt="n",xlab="",ylab="Species Richness")
legend(cex=3,"topright",legend=c("Wetlands","Controls"),fill=c("#D81B60","#1E88E5"))

### Coverage-based Shannon

covewetH <-abbcove[abbcove$Order.q==1,]$qD[c(1,3,5,7,9,11,13,15,17,19,21,23)]
covecontrH <-abbcove[abbcove$Order.q==1,]$qD[c(2,4,6,8,10,12,14,16,18,20,22,24)]
shapiro.test(covewetH-covecontrH)
t.test(covewetH,covecontrH,paired=T)
boxplot(covewetH,covecontrH)
coveH <-abbcove[abbcove$Order.q==1,]$qD
cove05H <-abbcove[abbcove$Order.q==1,]$qD.LCL
cove95H <-abbcove[abbcove$Order.q==1,]$qD.UCL
par(mar = c(2, 5, 2, 2) + 0.1)
plotCI(x=c(1:24),y=coveH,li=cove05H,ui=cove95H,lwd=2,cex.lab=2,cex.axis=1.5,pch=21,scol=c("#D81B60","#1E88E5"),xaxt="n",xlab="",ylab="Shannon Diversity")
legend(cex=3,"topright",legend=c("Wetlands","Controls"),fill=c("#D81B60","#1E88E5"))

### Coverage-based Simpson

covewetD <-abbcove[abbcove$Order.q==2,]$qD[c(1,3,5,7,9,11,13,15,17,19,21,23)]
covecontrD <-abbcove[abbcove$Order.q==2,]$qD[c(2,4,6,8,10,12,14,16,18,20,22,24)]
shapiro.test(covewetD-covecontrD)
t.test(covewetD,covecontrD,paired=T)
coveD <-abbcove[abbcove$Order.q==2,]$qD
cove05D <-abbcove[abbcove$Order.q==2,]$qD.LCL
cove95D <-abbcove[abbcove$Order.q==2,]$qD.UCL
par(mar = c(2, 5, 2, 2) + 0.1)
plotCI(x=c(1:24),y=coveD,li=cove05D,ui=cove95D,lwd=2,cex.lab=2,cex.axis=1.5,pch=21,scol=c("#D81B60","#1E88E5"),xaxt="n",xlab="",ylab="Simpson Diversity")
legend(cex=3,"topright",legend=c("Wetlands","Controls"),fill=c("#D81B60","#1E88E5"))


##############################TEST - FWS VS SEMI-NATURAL##############################


### Coverage-based Richness

covefwsS <-abbcove[abbcove$Order.q==0,]$qD[c(13,15,17,19,21,23)]
covenatS <-abbcove[abbcove$Order.q==0,]$qD[c(1,3,5,7,9,11)]
wilcox.test(covefwsS,covenatS)
boxplot(covefwsS,covenatS)
wetS <-abbcove[abbcove$Order.q==0,]$qD[c(1,3,5,7,9,11,13,15,17,19,21,23)]
wet05S <-abbcove[abbcove$Order.q==0,]$qD.LCL[c(1,3,5,7,9,11,13,15,17,19,21,23)]
wet95S <-abbcove[abbcove$Order.q==0,]$qD.UCL[c(1,3,5,7,9,11,13,15,17,19,21,23)]
par(mar = c(2, 5, 2, 2) + 0.1)
plotCI(x=c(1:12),y=wetS,li=wet05S,ui=wet95S,lwd=2,cex.lab=2,cex.axis=1.5,pch=21,scol=c(rep("#FFC107",6),rep("#004D40",6)),xaxt="n",xlab="",ylab="Species Richness")
legend(cex=3,"topright",legend=c("Semi-Natural","FWS"),fill=c("#FFC107","#004D40"))

### Coverage-based Shannon

covefwsH <-abbcove[abbcove$Order.q==1,]$qD[c(13,15,17,19,21,23)]
covenatH <-abbcove[abbcove$Order.q==1,]$qD[c(1,3,5,7,9,11)]
shapiro.test(covefwsH)
shapiro.test(covenatH)
var.test(covefwsH,covenatH)
t.test(covefwsH,covenatH,var=T)
boxplot(covefwsH,covenatH)
wetH <-abbcove[abbcove$Order.q==1,]$qD[c(1,3,5,7,9,11,13,15,17,19,21,23)]
wet05H <-abbcove[abbcove$Order.q==1,]$qD.LCL[c(1,3,5,7,9,11,13,15,17,19,21,23)]
wet95H <-abbcove[abbcove$Order.q==1,]$qD.UCL[c(1,3,5,7,9,11,13,15,17,19,21,23)]
par(mar = c(2, 5, 2, 2) + 0.1)
plotCI(x=c(1:12),y=wetH,li=wet05H,ui=wet95H,lwd=2,cex.lab=2,cex.axis=1.5,pch=21,scol=c(rep("#FFC107",6),rep("#004D40",6)),xaxt="n",xlab="",ylab="Shannon Diversity")
legend(cex=3,"topright",legend=c("Semi-Natural","FWS"),fill=c("#FFC107","#004D40"))

### Coverage-based Simpson

covefwsD <-abbcove[abbcove$Order.q==2,]$qD[c(13,15,17,19,21,23)]
covenatD <-abbcove[abbcove$Order.q==2,]$qD[c(1,3,5,7,9,11)]
shapiro.test(covefwsD)
shapiro.test(covenatD)
var.test(covefwsD,covenatD)
t.test(covefwsD,covenatD,var=T)
boxplot(covefwsD,covenatD)
wetD <-abbcove[abbcove$Order.q==2,]$qD[c(1,3,5,7,9,11,13,15,17,19,21,23)]
wet05D <-abbcove[abbcove$Order.q==2,]$qD.LCL[c(1,3,5,7,9,11,13,15,17,19,21,23)]
wet95D <-abbcove[abbcove$Order.q==2,]$qD.UCL[c(1,3,5,7,9,11,13,15,17,19,21,23)]
par(mar = c(2, 5, 2, 2) + 0.1)
plotCI(x=c(1:12),y=wetD,li=wet05D,ui=wet95D,lwd=2,cex.lab=2,cex.axis=1.5,pch=21,scol=c(rep("#FFC107",6),rep("#004D40",6)),xaxt="n",xlab="",ylab="Simpson Diversity")
legend(cex=3,"topright",legend=c("Semi-Natural","FWS"),fill=c("#FFC107","#004D40"))


##############################GLM FOR S WITH ENV+CLC##############################


### Create a dataframe for analysis

coveS
Sanalysis <-cbind(coveS,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,clc$hydro,clc$phrg)
Sanalysis <-as.data.frame(Sanalysis)
colnames(Sanalysis)<-c("Scove","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
Sanalysis
str(Sanalysis)

### Create a model for Scove

cor(Sanalysis) #NDVI seems to be correlated with Imperviousness and wood
modelScove1<-glm(Scove~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=Sanalysis,family=poisson)
car::vif(modelScove1) #wood, not NDVI, is the one with higher VIF
modelScove1<-glm(Scove~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=Sanalysis,family=poisson)
car::vif(modelScove1) #VIF<5 it's ok
summary(modelScove1)
### You shouldn't use a Poisson distribution for non-integer response values.
### The AIC is based on the negative log-likelihood, which in turn is based on the log probability of 
### the observed values given the model. The probability of a non-integer value is zero, so the 
### log-likelihood is -Inf, so the negative log-likelihood is Inf. If you have even a single 
### non-integer value in your data set, this will happen. 
### We can't use AIC like this...
modelScove2<-update(modelScove1,~. - phrg) 
anova(modelScove1,modelScove2,test="Chi")
summary(modelScove2)
modelScove3<-update(modelScove2,~. - Hclc)
anova(modelScove2,modelScove3,test="Chi")
summary(modelScove3)
modelScove4<-update(modelScove3,~. - NDVI)
anova(modelScove3,modelScove4,test="Chi")
summary(modelScove4)
modelScove5<-update(modelScove4,~. - agro)
anova(modelScove4,modelScove5,test="Chi")
summary(modelScove5)
modelScove6<-update(modelScove5,~. - hydro)
anova(modelScove5,modelScove6,test="Chi")
summary(modelScove6)

### Check assumptions for Scove

check_model(modelScove6)
# Linear relationship between the predictors and the log of the response variable
Scoveres<-residuals(modelScove6) 
Scoveyfit<-fitted(modelScove6)
plot(Scoveyfit,Scoveres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(Scoveres~Scoveyfit))
abline(h = 0, col = "red")
# Data should not exibit overdispersion, mean of y should be equal to its variance
dispersiontest(modelScove6)

(modelScove6$null.deviance-modelScove6$deviance)/modelScove6$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(Scove~Imperviousness,data=Sanalysis,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Species Richness",xlab="Imperviousness%")
abline(glm(Scove~Imperviousness,data=Sanalysis))


##############################GLM FOR H WITH ENV+CLC##############################


### Create a dataframe for analysis

coveH
Hanalysis <-cbind(coveH,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,clc$hydro,clc$phrg)
Hanalysis <-as.data.frame(Hanalysis)
colnames(Hanalysis)<-c("Hcove","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
Hanalysis
str(Hanalysis)

### Create a model for Hcove

cor(Hanalysis) #NDVI seems to be correlated with Imperviousness and wood
modelHcove1<-glm(Hcove~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=Hanalysis,family=gaussian)
car::vif(modelHcove1) #wood, not NDVI, is the one with higher VIF
modelHcove1<-glm(Hcove~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=Hanalysis,family=gaussian)
car::vif(modelHcove1) #VIF<5 it's ok
summary(modelHcove1)
modelHcove2<-step(modelHcove1)
anova(modelHcove1,modelHcove2,test="Chi")
summary(modelHcove2)

### Check assumptions for Hcove

check_model(modelHcove2)
# Linear relationship between the predictors and the response variable
Hcoveres<-residuals(modelHcove2) 
Hcoveyfit<-fitted(modelHcove2)
plot(Hcoveyfit,Hcoveres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(Hcoveres~Hcoveyfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(Hcoveres))
abline(0,1)
shapiro.test(Hcoveres) 
# Homoscedasticity of residuals
plot(Hcoveyfit,abs(Hcoveres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gHcove<-lm(abs(Hcoveres)~Hcoveyfit)
abline(gHcove)
summary(gHcove)

(modelHcove2$null.deviance-modelHcove2$deviance)/modelHcove2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(Hcove~Imperviousness,data=Hanalysis,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Shannon Diversity",xlab="Imperviousness%")
abline(glm(Hcove~Imperviousness,data=Hanalysis))


##############################GLM FOR D WITH ENV+CLC##############################


### Create a dataframe for analysis

coveD
Danalysis <-cbind(coveD,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,clc$hydro,clc$phrg)
Danalysis <-as.data.frame(Danalysis)
colnames(Danalysis)<-c("Dcove","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
Danalysis
str(Danalysis)

### Create a model for Dcove

cor(Danalysis) #NDVI seems to be correlated with Imperviousness and wood
modelDcove1<-glm(Dcove~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=Danalysis,family=gaussian)
car::vif(modelDcove1) #wood, not NDVI, is the one with higher VIF
modelDcove1<-glm(Dcove~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=Danalysis,family=gaussian)
car::vif(modelDcove1) #VIF<5 it's ok
summary(modelDcove1)
modelDcove2<-step(modelDcove1)
anova(modelDcove1,modelDcove2,test="Chi")
summary(modelDcove2)

### Check assumptions for Dcove

check_model(modelDcove2)
# Linear relationship between the predictors and the response variable
Dcoveres<-residuals(modelDcove2) 
Dcoveyfit<-fitted(modelDcove2)
plot(Dcoveyfit,Dcoveres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(Dcoveres~Dcoveyfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(Dcoveres))
abline(0,1)
shapiro.test(Dcoveres) 
# Homoscedasticity of residuals
plot(Dcoveyfit,abs(Dcoveres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gDcove<-lm(abs(Dcoveres)~Dcoveyfit)
abline(gDcove)
summary(gDcove)

(modelDcove2$null.deviance-modelDcove2$deviance)/modelDcove2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(Dcove~Imperviousness,data=Danalysis,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Simpson Diversity",xlab="Imperviousness%")
abline(glm(Dcove~Imperviousness,data=Danalysis))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(Dcove~phrg,data=Danalysis,col="#1E88E5",cex.axis=1.5,cex.lab=2,cex=1.5,pch=17,ylab="Simpson Diversity",xlab="Reedbed Coverage%")
abline(glm(Dcove~phrg,data=Danalysis))


##############################GUILDS INTRO##############################


### Load dataset

### Guilds
WAT <-read.table("Waterbirds.txt")
PAS <-read.table("Reedbed-Passerines.txt")
ERB <-read.table("Herbivores-Waders.txt")
CAR <-read.table("Carnivores.txt")

### IUCN Group
IUCN <-read.table("IUCN.txt")


### We can't use a R/E because there are some plots with very poor abundance of wetland birds
### Differences are caused by the presence of a water body, using a R/E we would lose this difference


##############################T-TEST GUILD WAT##############################


### We can compute diversity indices on the abundance data for WAT guild
Category <-c(3,4,3,4,3,4,3,4,3,4,3,4,1,2,1,2,1,2,1,2,1,2,1,2)
Category <-factor(Category)
SWAT <-hill_taxa(t(WAT), q = 0)
HWAT <-hill_taxa(t(WAT), q = 1)
DWAT <-hill_taxa(t(WAT), q = 2)
WATdiv <-cbind(Category,SWAT,HWAT,DWAT)
WATdiv <-as.data.frame(WATdiv)
WATdiv$Category <-factor(WATdiv$Category)
WATwet <-rbind(WATdiv[WATdiv$Category==1,],WATdiv[WATdiv$Category==3,])
WATcontr <-rbind(WATdiv[WATdiv$Category==2,],WATdiv[WATdiv$Category==4,])
WATwetnat <-WATdiv[WATdiv$Category==3,]
WATcontrnat <-WATdiv[WATdiv$Category==4,]
WATwetfws <-WATdiv[WATdiv$Category==1,]
WATcontrfws <-WATdiv[WATdiv$Category==2,]


### We make a t-test between wetlands and control plots
wilcox.test(WATwet$SWAT,WATcontr$SWAT,paired=T)              
boxplot(WATwet$SWAT,WATcontr$SWAT)

shapiro.test(WATwet$HWAT)
shapiro.test(WATcontr$HWAT)
t.test(WATwet$HWAT,WATcontr$HWAT,paired=T) #we can't use a var.test on paired data
boxplot(WATwet$HWAT,WATcontr$HWAT)

shapiro.test(WATwet$DWAT)
shapiro.test(WATcontr$DWAT)
t.test(WATwet$DWAT,WATcontr$DWAT,paired=T)
boxplot(WATwet$DWAT,WATcontr$DWAT)


### We make a t-test between semi-natural and FWS plots
wilcox.test(WATwetnat$SWAT,WATwetfws$SWAT)
boxplot(WATwetnat$SWAT,WATwetfws$SWAT)

shapiro.test(WATwetnat$HWAT)
shapiro.test(WATwetfws$HWAT)
var.test(WATwetnat$HWAT,WATwetfws$HWAT)
t.test(WATwetnat$HWAT,WATwetfws$HWAT,var.equal=T)
boxplot(WATwetnat$HWAT,WATwetfws$HWAT)

shapiro.test(WATwetnat$DWAT)
shapiro.test(WATwetfws$DWAT)
var.test(WATwetnat$DWAT,WATwetfws$DWAT)
t.test(WATwetnat$DWAT,WATwetfws$DWAT,var.equal=T)
boxplot(WATwetnat$DWAT,WATwetfws$DWAT)


par(mar = c(2, 5, 2, 2) + 0.1)
plot(WATdiv$SWAT~WATdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Species Richness")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(WATdiv$HWAT~WATdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Shannon Diversity")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(WATdiv$DWAT~WATdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Simpson Diversity")
legend(cex=2,"topright",legend=c("FWS","FWS Control","Semi-Natural","Semi-Natural Control"),fill=c("#1E88E5","#FFC107","#004D40","#D81B60"))


##############################T-TEST GUILD PAS##############################


### We can compute diversity indices on the abundance data for PAS guild
SPAS <-hill_taxa(t(PAS), q = 0)
HPAS <-hill_taxa(t(PAS), q = 1)
DPAS <-hill_taxa(t(PAS), q = 2)
PASdiv <-cbind(Category,SPAS,HPAS,DPAS)
PASdiv <-as.data.frame(PASdiv)
PASdiv$Category <-factor(PASdiv$Category)
PASwet <-rbind(PASdiv[PASdiv$Category==1,],PASdiv[PASdiv$Category==3,])
PAScontr <-rbind(PASdiv[PASdiv$Category==2,],PASdiv[PASdiv$Category==4,])
PASwetnat <-PASdiv[PASdiv$Category==3,]
PAScontrnat <-PASdiv[PASdiv$Category==4,]
PASwetfws <-PASdiv[PASdiv$Category==1,]
PAScontrfws <-PASdiv[PASdiv$Category==2,]


### We make a t-test between wetlands and control plots
wilcox.test(PASwet$SPAS,PAScontr$SPAS,paired=T)
boxplot(PASwet$SPAS,PAScontr$SPAS)

shapiro.test(PASwet$HPAS) #not normal
wilcox.test(PASwet$HPAS,PAScontr$HPAS,paired=T)
boxplot(PASwet$HPAS,PAScontr$HPAS)

shapiro.test(PASwet$DPAS) #not normal
wilcox.test(PASwet$DPAS,PAScontr$DPAS,paired=T)
boxplot(PASwet$DPAS,PAScontr$DPAS)

### We make a t-test between semi-natural and FWS plots
shapiro.test(PASwetnat$DPAS) #not normal
wilcox.test(PASwetnat$DPAS,PASwetfws$DPAS)
boxplot(PASwetnat$DPAS,PASwetfws$DPAS)

par(mar = c(2, 5, 2, 2) + 0.1)
plot(PASdiv$SPAS~PASdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Species Richness")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(PASdiv$HPAS~PASdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Shannon Diversity")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(PASdiv$DPAS~PASdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Simpson Diversity")
legend(cex=2,"topright",legend=c("FWS","FWS Control","Semi-Natural","Semi-Natural Control"),fill=c("#1E88E5","#FFC107","#004D40","#D81B60"))


##############################T-TEST GUILD ERB##############################


### We can compute diversity indices on the abundance data for ERB guild
SERB <-hill_taxa(t(ERB), q = 0)
HERB <-hill_taxa(t(ERB), q = 1)
DERB <-hill_taxa(t(ERB), q = 2)
ERBdiv <-cbind(Category,SERB,HERB,DERB)
ERBdiv <-as.data.frame(ERBdiv)
ERBdiv$Category <-factor(ERBdiv$Category)
ERBwet <-rbind(ERBdiv[ERBdiv$Category==1,],ERBdiv[ERBdiv$Category==3,])
ERBcontr <-rbind(ERBdiv[ERBdiv$Category==2,],ERBdiv[ERBdiv$Category==4,])
ERBwetnat <-ERBdiv[ERBdiv$Category==3,]
ERBcontrnat <-ERBdiv[ERBdiv$Category==4,]
ERBwetfws <-ERBdiv[ERBdiv$Category==1,]
ERBcontrfws <-ERBdiv[ERBdiv$Category==2,]


### We make a t-test between wetlands and control plots
wilcox.test(ERBwet$SERB,ERBcontr$SERB,paired=T)
boxplot(ERBwet$SERB,ERBcontr$SERB)

shapiro.test(ERBwet$HERB)
shapiro.test(ERBcontr$HERB) #not normal
wilcox.test(ERBwet$HERB,ERBcontr$HERB,paired=T)
boxplot(ERBwet$HERB,ERBcontr$HERB)

shapiro.test(ERBwet$DERB) #not normal
wilcox.test(ERBwet$DERB,ERBcontr$DERB,paired=T)
boxplot(ERBwet$DERB,ERBcontr$DERB)


### We make a t-test between semi-natural and FWS plots
wilcox.test(ERBwetnat$SERB,ERBwetfws$SERB)
boxplot(ERBwetnat$SERB,ERBwetfws$SERB)

shapiro.test(ERBwetnat$HERB)
shapiro.test(ERBwetfws$HERB)
var.test(ERBwetnat$HERB,ERBwetfws$HERB)
t.test(ERBwetnat$HERB,ERBwetfws$HERB,var.equal=T)
boxplot(ERBwetnat$HERB,ERBwetfws$HERB)

shapiro.test(ERBwetnat$DERB)
shapiro.test(ERBwetfws$DERB) #not normal
wilcox.test(ERBwetnat$DERB,ERBwetfws$DERB)
boxplot(ERBwetnat$DERB,ERBwetfws$DERB)


par(mar = c(2, 5, 2, 2) + 0.1)
plot(ERBdiv$SERB~ERBdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Species Richness")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(ERBdiv$HERB~ERBdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Shannon Diversity")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(ERBdiv$DERB~ERBdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Simpson Diversity")
legend(cex=2,"topright",legend=c("FWS","FWS Control","Semi-Natural","Semi-Natural Control"),fill=c("#1E88E5","#FFC107","#004D40","#D81B60"))


##############################T-TEST GUILD CAR##############################


### We can compute diversity indices on the abundance data for CAR guild
SCAR <-hill_taxa(t(CAR), q = 0)
HCAR <-hill_taxa(t(CAR), q = 1)
DCAR <-hill_taxa(t(CAR), q = 2)
CARdiv <-cbind(Category,SCAR,HCAR,DCAR)
CARdiv <-as.data.frame(CARdiv)
CARdiv$Category <-factor(CARdiv$Category)
CARwet <-rbind(CARdiv[CARdiv$Category==1,],CARdiv[CARdiv$Category==3,])
CARcontr <-rbind(CARdiv[CARdiv$Category==2,],CARdiv[CARdiv$Category==4,])
CARwetnat <-CARdiv[CARdiv$Category==3,]
CARcontrnat <-CARdiv[CARdiv$Category==4,]
CARwetfws <-CARdiv[CARdiv$Category==1,]
CARcontrfws <-CARdiv[CARdiv$Category==2,]


### We make a t-test between wetlands and control plots
wilcox.test(CARwet$SCAR,CARcontr$SCAR,paired=T)
boxplot(CARwet$SCAR,CARcontr$SCAR)

shapiro.test(CARwet$HCAR) #not normal
wilcox.test(CARwet$HCAR,CARcontr$HCAR,paired=T)
boxplot(CARwet$HCAR,CARcontr$HCAR)

shapiro.test(CARwet$DCAR) #not normal
wilcox.test(CARwet$DCAR,CARcontr$DCAR,paired=T)
boxplot(CARwet$DCAR,CARcontr$DCAR)


### We make a t-test between semi-natural and FWS plots
wilcox.test(CARwetnat$SCAR,CARwetfws$SCAR)
boxplot(CARwetnat$SCAR,CARwetfws$SCAR)

shapiro.test(CARwetnat$HCAR)
shapiro.test(CARwetfws$HCAR)
var.test(CARwetnat$HCAR,CARwetfws$HCAR)
t.test(CARwetnat$HCAR,CARwetfws$HCAR,var.equal=T)
boxplot(CARwetnat$HCAR,CARwetfws$HCAR)

shapiro.test(CARwetnat$DCAR)
shapiro.test(CARwetfws$DCAR) #not normal
wilcox.test(CARwetnat$DCAR,CARwetfws$DCAR)
boxplot(CARwetnat$DCAR,CARwetfws$DCAR)


par(mar = c(2, 5, 2, 2) + 0.1)
plot(CARdiv$SCAR~CARdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Species Richness")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(CARdiv$HCAR~CARdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Shannon Diversity")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(CARdiv$DCAR~CARdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Simpson Diversity")
legend(cex=2,"topright",legend=c("FWS","FWS Control","Semi-Natural","Semi-Natural Control"),fill=c("#1E88E5","#FFC107","#004D40","#D81B60"))


##############################T-TEST IUCN##############################


### We can compute diversity indices on the abundance data for IUCN group
SIUCN <-hill_taxa(t(IUCN), q = 0)
HIUCN <-hill_taxa(t(IUCN), q = 1)
DIUCN <-hill_taxa(t(IUCN), q = 2)
IUCNdiv <-cbind(Category,SIUCN,HIUCN,DIUCN)
IUCNdiv <-as.data.frame(IUCNdiv)
IUCNdiv$Category <-factor(IUCNdiv$Category)
IUCNwet <-rbind(IUCNdiv[IUCNdiv$Category==1,],IUCNdiv[IUCNdiv$Category==3,])
IUCNcontr <-rbind(IUCNdiv[IUCNdiv$Category==2,],IUCNdiv[IUCNdiv$Category==4,])
IUCNwetnat <-IUCNdiv[IUCNdiv$Category==3,]
IUCNcontrnat <-IUCNdiv[IUCNdiv$Category==4,]
IUCNwetfws <-IUCNdiv[IUCNdiv$Category==1,]
IUCNcontrfws <-IUCNdiv[IUCNdiv$Category==2,]


### We make a t-test between wetlands and control plots
wilcox.test(IUCNwet$SIUCN,IUCNcontr$SIUCN,paired=T)
boxplot(IUCNwet$SIUCN,IUCNcontr$SIUCN)

shapiro.test(IUCNwet$HIUCN) #not normal
wilcox.test(IUCNwet$HIUCN,IUCNcontr$HIUCN,paired=T)
boxplot(IUCNwet$HIUCN,IUCNcontr$HIUCN)

shapiro.test(IUCNwet$DIUCN) #not normal
wilcox.test(IUCNwet$DIUCN,IUCNcontr$DIUCN,paired=T)
boxplot(IUCNwet$DIUCN,IUCNcontr$DIUCN)


par(mar = c(2, 5, 2, 2) + 0.1)
plot(IUCNdiv$SIUCN~IUCNdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Species Richness")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(IUCNdiv$HIUCN~IUCNdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Shannon Diversity")
par(mar = c(2, 5, 2, 2) + 0.1)
plot(IUCNdiv$DIUCN~IUCNdiv$Category,cex.axis=1.5,cex.lab=2,col=c("#1E88E5","#FFC107","#004D40","#D81B60"),xaxt="n",xlab="",ylab="Simpson Diversity")
legend(cex=2,"topright",legend=c("FWS","FWS Control","Semi-Natural","Semi-Natural Control"),fill=c("#1E88E5","#FFC107","#004D40","#D81B60"))


##############################GLM FOR WAT WITH ENV+CLC##############################


### Create a dataframe for analysis

WATanalysisS <-cbind(SWAT,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
WATanalysisH <-cbind(HWAT,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
WATanalysisD <-cbind(DWAT,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
WATanalysisS <-as.data.frame(WATanalysisS)
WATanalysisH <-as.data.frame(WATanalysisH)
WATanalysisD <-as.data.frame(WATanalysisD)
colnames(WATanalysisS)<-c("S","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(WATanalysisH)<-c("H","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(WATanalysisD)<-c("D","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
WATanalysisS
WATanalysisH
WATanalysisD
str(WATanalysisS)
str(WATanalysisH)
str(WATanalysisD)

### Create a model for WAT

cor(WATanalysisS) #NDVI seems to be correlated with Imperviousness and wood
modelWATS1<-glm(S~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=WATanalysisS,family=poisson)
car::vif(modelWATS1) #wood, not NDVI, is the one with higher VIF
modelWATS1<-glm(S~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=WATanalysisS,family=poisson)
car::vif(modelWATS1) #VIF<5 it's ok
summary(modelWATS1)
modelWATS2<-step(modelWATS1)
anova(modelWATS1,modelWATS2,test="Chi")
summary(modelWATS2)

cor(WATanalysisH) #NDVI seems to be correlated with Imperviousness and wood
modelWATH1<-glm(H~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=WATanalysisH,family=gaussian)
car::vif(modelWATH1) #wood, not NDVI, is the one with higher VIF
modelWATH1<-glm(H~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=WATanalysisH,family=gaussian)
car::vif(modelWATH1) #VIF<5 it's ok
summary(modelWATH1)
modelWATH2<-step(modelWATH1)
anova(modelWATH1,modelWATH2,test="Chi")
summary(modelWATH2)

cor(WATanalysisD) #NDVI seems to be correlated with Imperviousness and wood
modelWATD1<-glm(D~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=WATanalysisD,family=gaussian)
car::vif(modelWATD1) #wood, not NDVI, is the one with higher VIF
modelWATD1<-glm(D~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=WATanalysisD,family=gaussian)
car::vif(modelWATD1) #VIF<5 it's ok
summary(modelWATD1)
modelWATD2<-step(modelWATD1)
anova(modelWATD1,modelWATD2,test="Chi")
summary(modelWATD2)

### Check assumptions for WAT

check_model(modelWATS2)
# Linear relationship between the predictors and the log of the response variable
SWATres<-residuals(modelWATS2) 
SWATyfit<-fitted(modelWATS2)
plot(SWATyfit,SWATres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(SWATres~SWATyfit))
abline(h = 0, col = "red")
# Data should not exibit overdispersion, mean of y should be equal to its variance
dispersiontest(modelWATS2)

(modelWATS2$null.deviance-modelWATS2$deviance)/modelWATS2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~Imperviousness,data=WATanalysisS,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Species Richness",xlab="Imperviousness%")
abline(glm(S~Imperviousness,data=WATanalysisS))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~hydro,data=WATanalysisS,col="#004D40",cex.axis=1.5,cex.lab=2,cex=1.5,pch=18,ylab="Species Richness",xlab="Waterbody Coverage%")
abline(glm(S~hydro,data=WATanalysisS))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~phrg,data=WATanalysisS,col="#1E88E5",cex.axis=1.5,cex.lab=2,cex=1.5,pch=17,ylab="Species Richness",xlab="Reedbed Coverage%")
abline(glm(S~phrg,data=WATanalysisS))

check_model(modelWATH2)
# Linear relationship between the predictors and the response variable
HWATres<-residuals(modelWATH2) 
HWATyfit<-fitted(modelWATH2)
plot(HWATyfit,HWATres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(HWATres~HWATyfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(HWATres))
abline(0,1)
shapiro.test(HWATres) 
# Homoscedasticity of residuals
plot(HWATyfit,abs(HWATres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gHWAT<-lm(abs(HWATres)~HWATyfit)
abline(gHWAT)
summary(gHWAT)

(modelWATH2$null.deviance-modelWATH2$deviance)/modelWATH2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(H~Imperviousness,data=WATanalysisH,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Shannon Diversity",xlab="Imperviousness%")
abline(glm(H~Imperviousness,data=WATanalysisH))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(H~phrg,data=WATanalysisH,col="#1E88E5",cex.axis=1.5,cex.lab=2,cex=1.5,pch=17,ylab="Shannon Diversity",xlab="Reedbed Coverage%")
abline(glm(H~phrg,data=WATanalysisH))

check_model(modelWATD2)
# Linear relationship between the predictors and the response variable
DWATres<-residuals(modelWATD2) 
DWATyfit<-fitted(modelWATD2)
plot(DWATyfit,DWATres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(DWATres~DWATyfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(DWATres))
abline(0,1)
shapiro.test(DWATres) 
# Homoscedasticity of residuals
plot(DWATyfit,abs(DWATres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gDWAT<-lm(abs(DWATres)~DWATyfit)
abline(gDWAT)
summary(gDWAT)

(modelWATD2$null.deviance-modelWATD2$deviance)/modelWATD2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(D~Imperviousness,data=WATanalysisD,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Simpson Diversity",xlab="Imperviousness%")
abline(glm(D~Imperviousness,data=WATanalysisD))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(D~phrg,data=WATanalysisD,col="#1E88E5",cex.axis=1.5,cex.lab=2,cex=1.5,pch=17,ylab="Simpson Diversity",xlab="Reedbed Coverage%")
abline(glm(D~phrg,data=WATanalysisD))


##############################GLM FOR PAS WITH ENV+CLC##############################


### Create a dataframe for analysis

PASanalysisS <-cbind(SPAS,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
PASanalysisH <-cbind(HPAS,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
PASanalysisD <-cbind(DPAS,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
PASanalysisS <-as.data.frame(PASanalysisS)
PASanalysisH <-as.data.frame(PASanalysisH)
PASanalysisD <-as.data.frame(PASanalysisD)
colnames(PASanalysisS)<-c("S","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(PASanalysisH)<-c("H","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(PASanalysisD)<-c("D","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
PASanalysisS
PASanalysisH
PASanalysisD
str(PASanalysisS)
str(PASanalysisH)
str(PASanalysisD)

### Create a model for PAS

cor(PASanalysisS) #NDVI seems to be correlated with Imperviousness and wood
modelPASS1<-glm(S~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=PASanalysisS,family=poisson)
car::vif(modelPASS1) #wood, not NDVI, is the one with higher VIF
modelPASS1<-glm(S~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=PASanalysisS,family=poisson)
car::vif(modelPASS1) #VIF<5 it's ok
summary(modelPASS1)
modelPASS2<-step(modelPASS1)
anova(modelPASS1,modelPASS2,test="Chi")
summary(modelPASS2)

cor(PASanalysisH) #NDVI seems to be correlated with Imperviousness and wood
modelPASH1<-glm(H~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=PASanalysisH,family=gaussian)
car::vif(modelPASH1) #wood, not NDVI, is the one with higher VIF
modelPASH1<-glm(H~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=PASanalysisH,family=gaussian)
car::vif(modelPASH1) #VIF<5 it's ok
summary(modelPASH1)
modelPASH2<-step(modelPASH1)
anova(modelPASH1,modelPASH2,test="Chi")
summary(modelPASH2)

#PASanalysisD<-PASanalysisD[c(-9,-19),] #to remove Inf values
cor(PASanalysisD) #NDVI seems to be correlated with Imperviousness and wood
modelPASD1<-glm(D~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=PASanalysisD,family=gaussian)
car::vif(modelPASD1) #wood, not NDVI, is the one with higher VIF
modelPASD1<-glm(D~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=PASanalysisD,family=gaussian)
car::vif(modelPASD1) #VIF<5 it's ok
summary(modelPASD1)
modelPASD2<-step(modelPASD1)
anova(modelPASD1,modelPASD2,test="Chi")
summary(modelPASD2)

### Check assumptions for PAS

check_model(modelPASH2)
# Linear relationship between the predictors and the response variable
HPASres<-residuals(modelPASH2) 
HPASyfit<-fitted(modelPASH2)
plot(HPASyfit,HPASres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(HPASres~HPASyfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(HPASres))
abline(0,1)
shapiro.test(HPASres) 
# Homoscedasticity of residuals
plot(HPASyfit,abs(HPASres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gHPAS<-lm(abs(HPASres)~HPASyfit)
abline(gHPAS)
summary(gHPAS)

(modelPASH2$null.deviance-modelPASH2$deviance)/modelPASH2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(H~phrg,data=PASanalysisH,col="#1E88E5",cex.axis=1.5,cex.lab=2,cex=1.5,pch=17,ylab="Shannon Diversity",xlab="Reedbed Coverage%")
abline(glm(H~phrg,data=PASanalysisH))


##############################GLM FOR ERB WITH ENV+CLC##############################


### Create a dataframe for analysis

ERBanalysisS <-cbind(SERB,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
ERBanalysisH <-cbind(HERB,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
ERBanalysisD <-cbind(DERB,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
ERBanalysisS <-as.data.frame(ERBanalysisS)
ERBanalysisH <-as.data.frame(ERBanalysisH)
ERBanalysisD <-as.data.frame(ERBanalysisD)
colnames(ERBanalysisS)<-c("S","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(ERBanalysisH)<-c("H","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(ERBanalysisD)<-c("D","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
ERBanalysisS
ERBanalysisH
ERBanalysisD
str(ERBanalysisS)
str(ERBanalysisH)
str(ERBanalysisD)

### Create a model for ERB

cor(ERBanalysisS) #NDVI seems to be correlated with Imperviousness and wood
modelERBS1<-glm(S~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=ERBanalysisS,family=poisson)
car::vif(modelERBS1) #wood, not NDVI, is the one with higher VIF
modelERBS1<-glm(S~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=ERBanalysisS,family=poisson)
car::vif(modelERBS1) #VIF<5 it's ok
summary(modelERBS1)
modelERBS2<-step(modelERBS1)
anova(modelERBS1,modelERBS2,test="Chi")
summary(modelERBS2)

cor(ERBanalysisH) #NDVI seems to be correlated with Imperviousness and wood
modelERBH1<-glm(H~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=ERBanalysisH,family=gaussian)
car::vif(modelERBH1) #wood, not NDVI, is the one with higher VIF
modelERBH1<-glm(H~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=ERBanalysisH,family=gaussian)
car::vif(modelERBH1) #VIF<5 it's ok
summary(modelERBH1)
modelERBH2<-step(modelERBH1)
anova(modelERBH1,modelERBH2,test="Chi")
summary(modelERBH2)

#ERBanalysisD<-ERBanalysisD[c(-4,-8,-10,-12,-16,-19,-20,-24),] #to remove Inf values
cor(ERBanalysisD) #NDVI seems to be correlated with Imperviousness and wood
modelERBD1<-glm(D~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=ERBanalysisD,family=gaussian)
car::vif(modelERBD1) #NDVI is the one with higher VIF
modelERBD1<-glm(D~Imperviousness+Hclc+wood+agro+hydro+phrg,data=ERBanalysisD,family=gaussian)
car::vif(modelERBD1) #agro is the one with higher VIF
modelERBD1<-glm(D~Imperviousness+Hclc+wood+hydro+phrg,data=ERBanalysisD,family=gaussian)
car::vif(modelERBD1) #VIF<5 it's ok
summary(modelERBD1)
modelERBD2<-step(modelERBD1)
anova(modelERBD1,modelERBD2,test="Chi")
summary(modelERBD2)

### Check assumptions for ERB

check_model(modelERBS2)
# Linear relationship between the predictors and the log of the response variable
SERBres<-residuals(modelERBS2) 
SERByfit<-fitted(modelERBS2)
plot(SERByfit,SERBres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(SERBres~SERByfit))
abline(h = 0, col = "red")
# Data should not exibit overdispersion, mean of y should be equal to its variance
dispersiontest(modelERBS2)

(modelERBS2$null.deviance-modelERBS2$deviance)/modelERBS2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~Imperviousness,data=ERBanalysisS,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Species Richness",xlab="Imperviousness%")
abline(glm(S~Imperviousness,data=ERBanalysisS))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~hydro,data=ERBanalysisS,col="#004D40",cex.axis=1.5,cex.lab=2,cex=1.5,pch=18,ylab="Species Richness",xlab="Waterbody Coverage%")
abline(glm(S~hydro,data=ERBanalysisS))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~phrg,data=ERBanalysisS,col="#1E88E5",cex.axis=1.5,cex.lab=2,cex=1.5,pch=17,ylab="Species Richness",xlab="Reedbed Coverage%")
abline(glm(S~phrg,data=ERBanalysisS))

check_model(modelERBH2)
# Linear relationship between the predictors and the response variable
HERBres<-residuals(modelERBH2) 
HERByfit<-fitted(modelERBH2)
plot(HERByfit,HERBres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(HERBres~HERByfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(HERBres))
abline(0,1)
shapiro.test(HERBres) 
# Homoscedasticity of residuals
plot(HERByfit,abs(HERBres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gHERB<-lm(abs(HERBres)~HERByfit)
abline(gHERB)
summary(gHERB)

(modelERBH2$null.deviance-modelERBH2$deviance)/modelERBH2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(H~Imperviousness,data=ERBanalysisH,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Shannon Diversity",xlab="Imperviousness%")
abline(glm(H~Imperviousness,data=ERBanalysisH))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(H~phrg,data=ERBanalysisH,col="#1E88E5",cex.axis=1.5,cex.lab=2,cex=1.5,pch=17,ylab="Shannon Diversity",xlab="Reedbed Coverage%")
abline(glm(H~phrg,data=ERBanalysisH))

check_model(modelERBD2)
# Linear relationship between the predictors and the response variable
DERBres<-residuals(modelERBD2) 
DERByfit<-fitted(modelERBD2)
plot(DERByfit,DERBres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(DERBres~DERByfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(DERBres))
abline(0,1)
shapiro.test(DERBres) 
# Homoscedasticity of residuals
plot(DERByfit,abs(DERBres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gDERB<-lm(abs(DERBres)~DERByfit)
abline(gDERB)
summary(gDERB)

(modelERBD2$null.deviance-modelERBD2$deviance)/modelERBD2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(D~Hclc,data=ERBanalysisD,col="#FFC107",cex.axis=1.5,cex.lab=2,cex=1.5,pch=15,ylab="Simpson Diversity",xlab="Landscape Patch Shannon Diversity")
abline(glm(D~Hclc,data=ERBanalysisD))


##############################GLM FOR CAR WITH ENV+CLC##############################


### Create a dataframe for analysis

CARanalysisS <-cbind(SCAR,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
CARanalysisH <-cbind(HCAR,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
CARanalysisD <-cbind(DCAR,env$Imperviousness_percentage,env$NDVI_mean,clc$Hclc,clc$wood,clc$agro,
             clc$hydro,clc$phrg)
CARanalysisS <-as.data.frame(CARanalysisS)
CARanalysisH <-as.data.frame(CARanalysisH)
CARanalysisD <-as.data.frame(CARanalysisD)
colnames(CARanalysisS)<-c("S","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(CARanalysisH)<-c("H","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
colnames(CARanalysisD)<-c("D","Imperviousness","NDVI","Hclc","wood","agro","hydro","phrg")
CARanalysisS
CARanalysisH
CARanalysisD
str(CARanalysisS)
str(CARanalysisH)
str(CARanalysisD)

### Create a model for CAR

cor(CARanalysisS) #NDVI seems to be correlated with Imperviousness and wood
modelCARS1<-glm(S~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=CARanalysisS,family=poisson)
car::vif(modelCARS1) #wood, not NDVI, is the one with higher VIF
modelCARS1<-glm(S~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=CARanalysisS,family=poisson)
car::vif(modelCARS1) #Imperviousness is the one with higher VIF
modelCARS1<-glm(S~NDVI+Hclc+agro+hydro+phrg,data=CARanalysisS,family=poisson)
car::vif(modelCARS1) #VIF<5 it's ok
summary(modelCARS1)
modelCARS2<-step(modelCARS1)
anova(modelCARS1,modelCARS2,test="Chi")
summary(modelCARS2)

cor(CARanalysisH) #NDVI seems to be correlated with Imperviousness and wood
modelCARH1<-glm(H~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=CARanalysisH,family=gaussian)
car::vif(modelCARH1) #wood, not NDVI, is the one with higher VIF
modelCARH1<-glm(H~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=CARanalysisH,family=gaussian)
car::vif(modelCARH1) #VIF<5 it's ok
summary(modelCARH1)
modelCARH2<-step(modelCARH1)
anova(modelCARH1,modelCARH2,test="Chi")
summary(modelCARH2)

#CARanalysisD<-CARanalysisD[c(-4,-8,-10,-12,-15,-16,-20,-22,-24),] #to remove Inf values
cor(CARanalysisD) #NDVI seems to be correlated with Imperviousness and wood
modelCARD1<-glm(D~Imperviousness+NDVI+Hclc+wood+agro+hydro+phrg,data=CARanalysisD,family=gaussian)
car::vif(modelCARD1) #wood, not NDVI, is the one with higher VIF
modelCARD1<-glm(D~Imperviousness+NDVI+Hclc+agro+hydro+phrg,data=CARanalysisD,family=gaussian)
car::vif(modelCARD1) #NDVI is the one with higher VIF
modelCARD1<-glm(D~Imperviousness+Hclc+agro+hydro+phrg,data=CARanalysisD,family=gaussian)
car::vif(modelCARD1) #VIF<5 it's ok
summary(modelCARD1)
modelCARD2<-step(modelCARD1)
anova(modelCARD1,modelCARD2,test="Chi")
summary(modelCARD2)

### Check assumptions for CAR

check_model(modelCARS2)
# Linear relationship between the predictors and the log of the response variable
SCARres<-residuals(modelCARS2) 
SCARyfit<-fitted(modelCARS2)
plot(SCARyfit,SCARres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(SCARres~SCARyfit))
abline(h = 0, col = "red")
# Data should not exibit overdispersion, mean of y should be equal to its variance
dispersiontest(modelCARS2)

(modelCARS2$null.deviance-modelCARS2$deviance)/modelCARS2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~NDVI,data=CARanalysisS,col="#FFC107",cex.axis=1.5,cex.lab=2,cex=1.5,pch=15,ylab="Species Richness",xlab="NDVI")
abline(glm(S~NDVI,data=CARanalysisS))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(S~hydro,data=CARanalysisS,col="#004D40",cex.axis=1.5,cex.lab=2,cex=1.5,pch=18,ylab="Species Richness",xlab="Waterbody Coverage%")
abline(glm(S~hydro,data=CARanalysisS))

check_model(modelCARH2)
# Linear relationship between the predictors and the response variable
HCARres<-residuals(modelCARH2) 
HCARyfit<-fitted(modelCARH2)
plot(HCARyfit,HCARres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(HCARres~HCARyfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(HCARres))
abline(0,1)
shapiro.test(HCARres) 
# Homoscedasticity of residuals
plot(HCARyfit,abs(HCARres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gHCAR<-lm(abs(HCARres)~HCARyfit)
abline(gHCAR)
summary(gHCAR) #VERY NOT-HOMOSCEDASTIC

check_model(modelCARD2)
# Linear relationship between the predictors and the response variable
DCARres<-residuals(modelCARD2) 
DCARyfit<-fitted(modelCARD2)
plot(DCARyfit,DCARres,ylab="Residuals",xlab="Fitted",main="Residuals vs fitted")
abline(lm(DCARres~DCARyfit)) 
abline(h = 0, col = "red")
# Normality of response variable and residuals
qqnorm(scale(DCARres))
abline(0,1)
shapiro.test(DCARres) 
# Homoscedasticity of residuals
plot(DCARyfit,abs(DCARres),ylab="Residuals",xlab="Fitted",main="Residuals in absolute value vs fitted")
gDCAR<-lm(abs(DCARres)~DCARyfit)
abline(gDCAR)
summary(gDCAR)

(modelCARD2$null.deviance-modelCARD2$deviance)/modelCARD2$null.deviance
par(mar = c(5, 5, 2, 2) + 0.1)
plot(D~Imperviousness,data=CARanalysisD,col="#D81B60",cex.axis=1.5,cex.lab=2,cex=1.5,pch=16,ylab="Simpson Diversity",xlab="Imperviousness%")
abline(glm(D~Imperviousness,data=CARanalysisD))
par(mar = c(5, 5, 2, 2) + 0.1)
plot(D~hydro,data=CARanalysisD,col="#004D40",cex.axis=1.5,cex.lab=2,cex=1.5,pch=18,ylab="Simpson Diversity",xlab="Waterbody Coverage%")
abline(glm(D~hydro,data=CARanalysisD))
