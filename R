
#setwd('/Users/colleenchen/randomForest/post-OHBM/')
library(randomForest)
library(cluster)
library(gdata)
library(shiny)
set.seed(420)
### SAVE the variables and workspace !!!!
save.image("/Volumes/jedi/Colleen/post-OHBM/Unsupervised/Unsupervised.RData")
#save.image("~/randomForest/post-OHBM/Unsupervised/Unsupervised.RData")

# load("/Volumes/jedi/Colleen/post-OHBM/Unsupervised/Unsupervised.RData")

load("~/R/Unsupervised.RData")

########### 11/10/15 #################################################################################
age.reg.srf = randomForest(x=regressed, y=response, data=regressed,importance=TRUE,ntree=25000)
# OOB age regressed : 35.19%

age.reg.urf = randomForest(x=regressed, data=regressed,proximity=TRUE, importance=TRUE,ntree=25000)
MDSplot(age.reg.urf,fac=response,pch=as.numeric(response))





########### 10/26/15 #################################################################################
ABIDE449HO <- read.table("~/R/ABIDE_n449_HO.txt", quote="\"", comment.char="")
demo449HO = read.csv("~/R/Data_anat.csv", header = TRUE, sep=",")

response = as.factor(demo449HO$Diagnostic)

abide449ho.srf = randomForest(x=ABIDE449HO, y=response, data=ABIDE449HO, importance=TRUE,ntree=50000)
# OOB: 36.53% 
plot(abide449ho.srf)


regressed <- mat.or.vec(dim(ABIDE449HO)[1], dim(ABIDE449HO)[2])
# rdata <-      mat.or.vec(dim(ABIDE449HO)[1], dim(ABIDE449HO)[2])
# REGRESSED OUT AGE!!! 
for( i in 1: dim(ABIDE449HO)[2] ){
  lm.r = lm( ABIDE449HO[,i] ~ demo449HO$AGE ,data=ABIDE449HO )
  row = resid(lm.r)
  col = t(row) 
  regressed[,i] = t(col)  
  print(i)
}
dim(regressed)

age.reg.srf = randomForest(x=regressed, y=response, data=regressed,importance=TRUE,ntree=25000)
# OOB age regressed : 35.19%


avg.abide = NULL
avg.age=NULL
nsubj = dim(ABIDE449HO)[1]
for (i in 1:nsubj){
  avg.abide[i] = mean(as.matrix(ABIDE449HO[i,]))
  avg.age[i]= mean(regressed[i,])
  # avg.site[i]=mean(regressed.site2[i,])
}
par(mfrow=c(1,1))
boxplot(avg.abide ~ demo449HO$SITE, data=demo449HO)
boxplot(avg.age ~ demo449HO$SITE, data=demo449HO)
# boxplot(avg.site ~ demo$SITE, data=demo)



############################################################################################
### LOAD the new data after adjusting for number of ASD/TD per SITE 
##### NEW SAMPLE N=449
abide449 = read.table("ABIDE_n449.txt",header=FALSE,sep=" ")
##### NEW SAMPLE N=449HO
abide449ho = read.table("/Volumes/jedi/Colleen/HarvardOxfordROIs/ABIDE_n449_HO.txt",header=FALSE,sep=" ")
demo449ho = read.xls("/Volumes/jedi/Colleen/HarvardOxfordROIs/Data_anatomical1.xlsx", sheet = 1, header = TRUE)
# read in the deamographic information
demo = read.xls("DATA_site.xlsx", sheet = 1, header = TRUE)

abide449.srf = randomForest(x=abide449, y=as.factor(demo$Diagnostic), data=abide449, importance=TRUE,ntree=10000)
# oob: 39.4%

########### 9/25/15 #################################################################################

sites = factor(demo$SITE)
regressed <- mat.or.vec(dim(abide449)[1], dim(abide449)[2])
rdata <- mat.or.vec(dim(abide449)[1], dim(abide449)[2])
# REGRESSED OUT AGE!!! 
for( i in 1: dim(abide449)[2] ){
  lm.r = lm( abide449[,i] ~ demo$AGE ,data=abide449 )
   row = resid(lm.r)
   col = t(row) 
  regressed[,i] = t(col)  
  print(i)
}
# then regressed out SITE in regression-copy.R script 
ageSite.rf = randomForest(x=regressed.site2, y=as.factor(demo$Diagnostic),data=regressed.site2,ntree=10000,importance=TRUE)
# OOB 38%
############################# OOB 35%
varImpPlot(ageSite.rf,type=1)
predict(ageSite.rf)
agesite.mda=importance(ageSite.rf,type=1)
agesite.mda=order(agesite.mda,decreasing=TRUE)
features = agesite.mda[1:200]

agesite.urf = randomForest(x=regressed.site2[,features], data=regressed.site2,ntree=10000,importance=TRUE,proximity=TRUE)
pamSilhouetteWidth.fun(1 - agesite.urf$proximity)
pamData <- pam(1 - agesite.urf$proximity, k = 2, diss = TRUE)
mds.abide = MDSplot(agesite.urf, pamData$clustering, k=2, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
legend("bottomright",c("cluster 1","cluster 2"),pch=c(21:25)[unclass=pamData$clustering],col = rep("black", 3), pt.bg = rainbow(3), cex=.4)
text(mds.abide$points[,1]+.001,mds.abide$points[,2]-.002,as.character(demo$Diagnostic),cex=.6)
title(main="Unsupervised RF using top 100 features, PAM Groupings, k=3")
# record the PAM clustering results: 
ifelse(as.matrix(pamData$clustering)==1,1,0)
########################### Unsupervised accuracy = 65% 

ASD = regressed.site2[demo$Diagnostic==1,]
TD = regressed.site2[demo$Diagnostic==0,]

# how many groups in TD? 
td.urf = randomForest(x=TD[,features],data=regressed.site2,ntree=10000,importance=TRUE,proximity=TRUE)
pamSilhouetteWidth.fun(1 - td.urf$proximity)
pamData <- pam(1 - td.urf$proximity, k = 2, diss = TRUE)
mds.abide = MDSplot(td.urf, pamData$clustering, k=2, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
legend("bottomright",c("cluster 1","cluster 2"),pch=c(21:25)[unclass=pamData$clustering],col = rep("black", 3), pt.bg = rainbow(3), cex=.4)
# record the PAM clustering results: 
# ifelse(as.matrix(pamData$clustering)==1,1,0)

asd.urf = randomForest(x=ASD[,features],data=regressed.site2,ntree=10000,importance=TRUE,proximity=TRUE)
pamSilhouetteWidth.fun(1 - asd.urf$proximity)
pamData <- pam(1 - asd.urf$proximity, k = 2, diss = TRUE)
mds.abide = MDSplot(td.urf, pamData$clustering, k=2, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
# this was aweful! 


sitevar.rf = randomForest(x=regressed,y=demo$SITE,data=regressed,ntree=10000,importance=TRUE)
# OOB 62.58%
varImpPlot(sitevar.rf,type=1)


demo.rf = randomForest(x=demo[,4:22],y=as.factor(demo$Diagnostic),data=demo,ntree=10000,importance=TRUE)
demo.urf = randomForest(x=demo[,4:22],data=demo,ntree=10000,importance=TRUE,proximity=TRUE)




########### Axel 9/17 #################################################################################
as.factor(demo$SITE)
site=table(demo$SITE,demo$Diagnostic)
chisq.test(x=site)
# The null hypothesis is that there is no association between SITE and Diagnostic group 
# reject null if p-value is < 0.05
fisher.test(table(demo$SITE=='TRINITY',demo$Diagnostic))
fisher.test(table(demo[demo$SITE=='YALE',1],demo[demo$SITE=='YALE',3]))


###########after meeting with Axel after BM #################################################################################
# rdata is the AGE and SITE regressed data: 
rdata = mat.or.vec(nsubj, number_columns)
rdata[,1] = ifelse(regressed[,1]==2,1,0) # this turns the labels into 'characters' 
rdata = regressed[,2:number_columns] # NO DIAGNOSTIC LABELS! use "response" variable 
dim(rdata)
# supervised learning and feature selection using REGRESSED data 
rsrf = randomForest(x=rdata,y=as.factor(demo$Diagnostic), data=rdata,importance=TRUE,ntree=10000)
rsrf
#OOB estimate of  error rate: 35.06%
#Confusion matrix:
#  0   1 class.error
#0 173  67   0.2791667
#1  95 127   0.4279279
plot(rsrf)
varImpPlot(rsrf,type=1,n=200)
rmda = importance(rsrf,type=1)
rfeatures = order(rmda,decreasing=TRUE)

n=200
trees=10000

rsrf.200 = randomForest(x=rdata[,features[1:n]],y=response,data=rdata,importance=TRUE,ntree=trees)
rsrf.200
# 200 features, OOB: 22.29% 

## unsupervised learning of abide 

urf.abide.200 = randomForest(x=rdata[,features[1:n]],data=rdata,importance=TRUE,proximity=TRUE,ntree=trees)
pamSilhouetteWidth.fun(1 - urf.abide.200$proximity)
pamData <- pam(1 - urf.abide.200$proximity, k = 2, diss = TRUE)
mds.abide = MDSplot(urf.abide.200, pamData$clustering, k=2, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
legend("bottomright",c("cluster 1","cluster 2"),pch=c(21:25)[unclass=pamData$clustering],col = rep("black", 3), pt.bg = rainbow(3), cex=.4)
text(mds.abide$points[,1]+.001,mds.abide$points[,2]-.002,as.character(demo$Diagnostic),cex=.6)
title(main="Unsupervised RF using top 100 features, PAM Groupings, k=3")
# record the PAM clustering results: 
as.matrix(pamData$clustering)

library(scatterplot3d)
scatterplot3d(mds.abide$points, highlight.3d=TRUE, col.axis="blue", col.grid="lightblue",pch=c(21:25)[unclass=pamData$clustering],bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering])
text(mds.abide$points[,1], mds.abide$points[,2], mds.abide$points[,3]+.004, as.character(demo$Diagnostic), cex=.7)




#########after meeting with Barb before BM ##############################################################################

abide.urf = randomForest(x=allsubj1,data=allsubj1,importance=TRUE,proximity=TRUE)
# varImpPlot() distinguishing between the real dataset 
# and the synthetic dataset (the noisy version of the real dataset)
varImpPlot(abide.urf,type=1,n=5000)
varImpPlot(abide.urf,type=1,n=700)
mda.urf = importance(abide.urf,type=1)
n1=700
features.urf = order(mda.urf,decreasing=TRUE)[1:n1]
pamSilhouetteWidth.fun(1 - abide.urf$proximity)
# k=2 groups 
pamData <- pam(1 - abide.urf$proximity, k = 2, diss = TRUE)
mds=MDSplot(abide.urf, k=2, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",2))
r0 =as.factor(as.matrix(pamData$clustering))


urf.700 = randomForest(x=allsubj1[,features.urf],data=allsubj1,proximity=TRUE,importance=TRUE,ntree=10000)
varImpPlot(urf.700,type=1,n=700)
mda.urf = importance(urf.700,type=1)
features.urf = order(mda.urf,decreasing=TRUE)[1:200]
pamSilhouetteWidth.fun(1 - urf.700$proximity)
pamData <- pam(1 - urf.700$proximity, k = 3, diss = TRUE)
mds=MDSplot(urf.700, k=2, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",2))
as.matrix(pamData$clustering)

srf.700 = randomForest(x=allsubj1[,features.urf],y=r0,data=allsubj1,proximity=TRUE,importance=TRUE,ntree=10000)
srf.700



urf.200 = randomForest(x=allsubj1[,features.urf],data=allsubj1,proximity=TRUE,importance=TRUE,ntree=10000)

pamSilhouetteWidth.fun(1 - urf.200$proximity)
pamData <- pam(1 - urf.200$proximity, k = 2, diss = TRUE)
mds=MDSplot(urf.200, k=2, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",2))
text(mds$points[,1]+.001,mds$points[,2]-.002,as.character(demo$Diagnostic),cex=.6)
as.matrix(pamData$clustering)

presp = as.factor(as.matrix(pamData$clustering))
srf.urf.200 = randomForest(x=allsubj1[features.urf],y=presp,data=allsubj1,proximity=TRUE,importance=TRUE,ntree=10000)
srf.urf.200
# OOB: 22.5%
varImpPlot(srf.urf.200,type=1)

pamSilhouetteWidth.fun(1 - srf.urf.200$proximity)
pamData <- pam(1 - srf.urf.200$proximity, k = 2, diss = TRUE)
mds=MDSplot(srf.urf.200, k=2, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",2))
text(mds$points[,1]+.001,mds$points[,2]-.002,as.character(demo$Diagnostic),cex=.6)
as.matrix(pamData$clustering)

 



legend("bottomright",c("Cluster 1", "Cluster 2"), pch=c(21:22), col = rep("black", 2), pt.bg = rainbow(2), cex=.6)
library(scatterplot3d)
scatterplot3d(mds$points, highlight.3d=TRUE, col.axis="blue", col.grid="lightblue",pch=c(21:25)[unclass=pamData$clustering],bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering])

# K=2
pamData <- pam(1 - urf.asd$proximity, k = 2, diss = TRUE)
mds <- MDSplot(urf.asd, k=3, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
library(scatterplot3d)
scatterplot3d(mds$points, highlight.3d=TRUE, col.axis="blue", col.grid="lightblue",pch=c(21:25)[unclass=pamData$clustering],bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering])
text(mds$points[,1], mds$points[,2]-.004, as.character(demo$SUB_ID), cex=.7)

# ASD only using the 200 features selected from unsupervised variable importance MDA 
asd.50 = randomForest(x=ASD[features.urf],data=ASD,proximity=TRUE,importance=TRUE,ntree=10000)
pamSilhouetteWidth.fun(1 - asd.50$proximity)
pamData <- pam(1 - asd.50$proximity, k = 3, diss = TRUE)
mds=MDSplot(asd.50, k=3, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",2))
legend("bottomright",c("Cluster 1", "Cluster 2"), pch=c(21:22), col = rep("black", 2), pt.bg = rainbow(2), cex=.6)
library(scatterplot3d)
scatterplot3d(mds$points, highlight.3d=TRUE, col.axis="blue", col.grid="lightblue",pch=c(21:25)[unclass=pamData$clustering],bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering])


############################################################################################
# Demographic information of the sample
sample <- read.csv("IncraseSampleSize.csv",header=TRUE,sep=",")
attach(sample)
# AGE: 
summary( sample[ sample$Diagnostic==0, 4] )
sd (sample[ sample$Diagnostic==0, 4])
summary( sample[ sample$Diagnostic==1, 4] )
sd (sample[ sample$Diagnostic==1, 4])
# MOTION 
summary( sample[ sample$Diagnostic==0, 5] )
sd (sample[ sample$Diagnostic==0, 5])
summary( sample[ sample$Diagnostic==1, 5] )
sd (sample[ sample$Diagnostic==1, 5])
# TIMEPOINTS 
summary( sample[ sample$Diagnostic==0, 6] )
sd (sample[ sample$Diagnostic==0, 6])
summary( sample[ sample$Diagnostic==1, 6] )
sd (sample[ sample$Diagnostic==1, 6])
# NVIQ 
summary( sample[ sample$Diagnostic==0, 11] )
sd (sample[ sample$Diagnostic==0, 11])
summary( sample[ sample$Diagnostic==1, 11] )
sd (sample[ sample$Diagnostic==1, 11])

sample[sample==-9999]<-NA
# ADOS_TOTAL 12
summary( sample[ sample$Diagnostic==1, 12] )
sd (sample[ sample$Diagnostic==1, 12])

# ================= read in data! 
#setwd('/Volumes/Turing/UnsupervisedABIDE/Unsupervised')
# read in the data matrix ABIDE_n462
allsubj = read.table("ABIDE_n462.txt",header=FALSE,sep=" ")
# read in the deamographic information
demo = read.xls("DATA.xlsx", sheet = 1, header = TRUE)
# read in the ROI information
ROI = read.xls("DATA.xlsx", sheet = 2, header = TRUE)

################################################################# SUPERVISED binary (ASD vs TD) classification 
set.seed(420)
response = as.factor(demo$Diagnostic)
nsubj = dim(allsubj)[1]
nfeature = dim(allsubj)[2]
srf = randomForest(x=allsubj1,y=as.factor(demo$Diagnostic), data=allsubj1,importance=TRUE,ntree=10000)
srf
plot(srf)
varImpPlot(srf,type=1,n=100)
mda = importance(srf,type=1)
features = order(mda,decreasing=TRUE)

n=200
trees=10000

srf.200 = randomForest(x=allsubj1[,features[1:n]],y=response,data=allsubj1,importance=TRUE,ntree=trees)
srf.200
# 100 features, OOB: 20.13%
# 200 features, OOB: 19.26% 

#########----------------Random Forest ROC curve for SUPERVISED Learning 
library(ROCR)

str(allsubj)

allsubj = cbind(response,allsubj[,])

set1<-allsubj[allsubj$response=="1",] # ASD
set2<-allsubj[allsubj$response=="0",] # TD


asd = dim(set1)[1]
td = dim(set2)[1]

t.asd = asd *2/3  # 222 ASD, 58311   t.asd=148
t.td = td *2/3  # 240 TD, 58311      t.td=160

set.seed(4)
training1<-sample(1:asd,t.asd)
test1<-(1:asd)[-training1]
sum((1:asd)==sort(c(training1,test1))) # 222

training2<-sample(1:td,t.td)
test2<-(1:td)[-training2]
sum((1:td==sort(c(training2,test2)))) # 240

train <- rbind(set1[training1,],set2[training2,])
test <- rbind(set1[test1,], set2[test2,])

dim(train)  # 308 (148+160)  58312
dim(test)   #  154 (462-308)  58312

set.seed(420)
# response variable in training data set: ASD/TD
responsetrain <- train[,1]

# predictors/features in training data set: correlations 
predstrain <- train[,-1]

ABIDE.rf.train <- randomForest(x=predstrain,y=responsetrain, data=train, mtry=100, ntree=10000, importance=TRUE)
mda.train = importance(ABIDE.rf.train,type=1)
features.train = order(mda.train,decreasing=TRUE)


ABIDE.rf.pr = predict(ABIDE.rf.train, test, type = "prob")[,2] 
ABIDE.rf.pred = prediction(ABIDE.rf.pr, test$response)
ABIDE.rf.perf = performance(ABIDE.rf.pred,"tpr","fpr")
ABIDE.auc <- performance(ABIDE.rf.pred,"auc")
ABIDE.auc <- unlist(slot(ABIDE.auc, "y.values"))

#BUILDING A ROC curve with 100 features
n=100
f = features.train[1:n]
ABIDE100.rf <-randomForest(x=predstrain[,f],y=responsetrain, data=train, ntree=10000, mtry=3, importance=TRUE)

ABIDE100.rf.pr = predict(ABIDE100.rf, test, type = "prob")[,2] 
ABIDE100.rf.pred = prediction(ABIDE100.rf.pr, test$response)
ABIDE100.rf.perf = performance(ABIDE100.rf.pred,"tpr","fpr")
ABIDE100.auc <- performance(ABIDE100.rf.pred,"auc")
ABIDE100.auc <- unlist(slot(ABIDE100.auc, "y.values"))

# using features from training set... 
abide100.rf.pr = predict(ABIDE100.rf,test,type="prob")[,2]
abide100.rf.pred = prediction(abide100.rf.pr, test$response)
abide100.auc = performance(abide100.rf.pred,"auc")
abide100.auc = unlist(slot(abide100.auc, "y.values"))


# AUC display
minauc<-min(round(ABIDE.auc, digits = 2))
maxauc<-max(round(ABIDE100.auc, digits = 2))
minauct <- paste(c("min(AUC)  = "),minauc,sep="")
maxauct <- paste(c("max(AUC) = "),maxauc,sep="")

## plotting the ROC curves and AUC 
plot(ABIDE.rf.perf,main="ROC Curve for supervised Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
par(new=TRUE)
plot(ABIDE100.rf.perf,col=3,lwd=2)
legend("bottomright",c("all features: AUC 0.64","top 100 features: AUC 0.88"), pch=c(23:23), col =c(2,3), pt.bg = rainbow(3), cex=.8)
legend("topleft",c(minauct,maxauct),cex=0.7,box.col = "white",border="black")
pdf("/Users/colleenchen/randomForest/post-OHBM/figures/srf_ROC.pdf", height=6, width=6)
#legend(0.3,0.6,c(minauct,maxauct,"\n"),border="white",cex=1.7,box.col = "white")
#test.classRF <- predict (ABIDE100.rf, newdata=test, type='prob')
par(new=FALSE)
#################-----------UNSUPERVISED cluster of ASD vs TD using top features ---------------------------------------------

allsubj1 = allsubj[,-1]
urf.abide = randomForest(x=allsubj1[,features[1:200]],data=allsubj1,proximity=TRUE,ntree=10000 )

pamSilhouetteWidth.fun(1 - urf.abide$proximity)
pamData <- pam(1 - urf.abide$proximity, k = 2, diss = TRUE)
mds.abide = MDSplot(urf.abide, pamData$clustering, k=2, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
legend("bottomright",c("cluster 1","cluster 2","cluster 3"),pch=c(21:25)[unclass=pamData$clustering],col = rep("black", 3), pt.bg = rainbow(3), cex=.4)
text(mds.abide$points[,1]+.001,mds.abide$points[,2]-.002,as.character(demo$Diagnostic),cex=.6)
title(main="Unsupervised RF using top 100 features, PAM Groupings, k=3")
# record the PAM clustering results: 
as.matrix(pamData$clustering)

library(scatterplot3d)
scatterplot3d(mds.abide$points, highlight.3d=TRUE, col.axis="blue", col.grid="lightblue",pch=c(21:25)[unclass=pamData$clustering],bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering])
text(mds.abide$points[,1], mds.abide$points[,2], mds.abide$points[,3]+.004, as.character(demo$Diagnostic), cex=.7)

### clustering of ASD only using supervised features100
urf.asd200 = randomForest(x=ASD[,features[1:200]],data=ASD,proximity=TRUE,ntree=10000)
pamSilhouetteWidth.fun(1 - urf.asd200$proximity)
pamData <- pam(1 - urf.asd200$proximity, k = 2, diss = TRUE)
MDSplot(urf.asd100, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
legend("bottomright",c("cluster 1","cluster 2","cluster 3"),pch=c(21:25)[unclass=pamData$clustering],col = rep("black", 3), pt.bg = rainbow(3), cex=.4)
text(mds.urf.asd100$points[,1]+.001,mds.urf.asd100$points[,2]-.002,as.character(demo$Diagnostic),cex=.4)
title(main="Unsupervised RF using top 100 features, PAM Groupings, k=3")
# record the PAM clustering results: 
as.matrix(pamData$clustering)


##################################### VARIANCE FILTER DIMENTIONALITY REDUCTION ::
variance = mat.or.vec(nfeature,1)
for (i in 1:nfeature) {
  variance[i] = var(allsubj[,i])  
}
summary(variance)
#Min.     1st Qu.  Median     Mean    3rd Qu.    Max. 
#0.01016  0.03663  0.04109    0.04292 0.04739    0.11060 
plot(variance)
threshold = 0.08

features_hivar = which(variance>threshold) 
srf_hivar = randomForest(x=allsubj[,features_hivar],y=as.factor(demo$Diagnostic), data=allsubj,importance=TRUE,ntree=1000)
plot(srf_hivar)

urf_hivar = randomForest(x=allsubj[,features_hivar],data=allsubj,proximity=TRUE,importance=TRUE,ntree=10000)
MDSplot(urf_hivar, as.factor(demo$Diagnostic),  pch=demo$Diagnostic)
pamSilhouetteWidth.fun(1 - urf_hivar$proximity)
pamData <- pam(1 - urf_hivar$proximity, k = 2, diss = TRUE)
MDSplot(urf_hivar, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))

#################################### split into ASD data 

ASD = allsubj[demo$Diagnostic==1,]
TD = allsubj[demo$Diagnostic==0,]

urf.td = randomForest(x=TD[,features[1:n]],data=TD,proximity=TRUE,importance=TRUE,ntree=1000)

pamSilhouetteWidth.fun(1 - urf.td$proximity)
pamData <- pam(1 - urf.td$proximity, k = 2, diss = TRUE)
MDSplot(urf.td, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
legend("bottomright",c("Cluster 1", "Cluster 2","Cluster 3"),pch=c(21:23), col = rep("black", 3), pt.bg = rainbow(3), cex=.8)


urf.asd = randomForest(x=ASD[,features_hivar],data=ASD,proximity=TRUE,importance=TRUE,ntree=1000)
pamSilhouetteWidth.fun(1 - urf.asd$proximity)
pamData <- pam(1 - urf.asd$proximity, k = 2, diss = TRUE)
MDSplot(urf.asd, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",2))
legend("bottomright",c("Cluster 1", "Cluster 2"), pch=c(21:22), col = rep("black", 2), pt.bg = rainbow(2), cex=.6)

# K=2
pamData <- pam(1 - urf.asd$proximity, k = 2, diss = TRUE)
mds <- MDSplot(urf.asd, k=3, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))

library(scatterplot3d)
scatterplot3d(mds$points, highlight.3d=TRUE, col.axis="blue", col.grid="lightblue",pch=c(21:25)[unclass=pamData$clustering],bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering])
text(mds$points[,1], mds$points[,2]-.004, as.character(demo$SUB_ID), cex=.7)
#!!!!!!!!!!!! Scatter plot 3D. a pams function call, which can give you k=3

#=== Lower variance threshold 
lo_thresh = 0.04
features_lovar = which(variance>lo_thresh) 
lo.urf.asd = randomForest(x=ASD[,features_lovar],data=ASD,proximity=TRUE,importance=TRUE,ntree=1000)

pamSilhouetteWidth.fun(1 - lo.urf.asd$proximity)
pamData <- pam(1 - lo.urf.asd$proximity, k = 2, diss = TRUE)
MDSplot(lo.urf.asd, pamData$clustering, bg=rainbow(length(levels(as.factor(pamData$clustering))))[unclass=pamData$clustering], pch=c(21:25)[unclass=pamData$clustering], palette=rep("black",5))
legend("bottomright",c("Cluster 1", "Cluster 2","Cluster 3"), 
       pch=c(21:23), col = rep("black", 3), pt.bg = rainbow(3), cex=.8)


### SAVE the variables and workspace !!!!
save.image("/Volumes/jedi/Colleen/post-OHBM/Unsupervised/Unsupervised.RData")
#save.image("~/randomForest/post-OHBM/Unsupervised/Unsupervised.RData")


 
