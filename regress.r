
#setwd('/Volumes/Purkinje/SVM-RFE/')
#source('regression.r')

#heisenberg <- read.csv(file="ABIDE_SDSU_demo.csv",head=TRUE,sep=",")
heisenberg = demo
nsubj=dim(demo)[1]
a<-rep(1, nsubj)
heisenberg[,"site_dumb01"] = a
heisenberg[,"site_dumb02"] = a
heisenberg[,"site_dumb03"] = a
heisenberg[,"site_dumb04"] = a
heisenberg[,"site_dumb05"] = a
heisenberg[,"site_dumb06"] = a
heisenberg[,"site_dumb07"] = a
heisenberg[,"site_dumb08"] = a
heisenberg[,"site_dumb09"] = a
heisenberg[,"site_dumb10"] = a
heisenberg[,"site_dumb11"] = a
heisenberg[,"site_dumb12"] = a
heisenberg[,"site_dumb13"] = a
heisenberg <- within(heisenberg, site_dumb01<-ifelse(SITE=="CMU", 1, 0))
heisenberg <- within(heisenberg, site_dumb02<-ifelse(SITE=="KKI", 1, 0))
heisenberg <- within(heisenberg, site_dumb03<-ifelse(SITE=="LEUVEN_1", 1, 0))
heisenberg <- within(heisenberg, site_dumb04<-ifelse(SITE=="LEUVEN_2", 1, 0))
heisenberg <- within(heisenberg, site_dumb05<-ifelse(SITE=="NYU", 1, 0))
heisenberg <- within(heisenberg, site_dumb06<-ifelse(SITE=="OLIN", 1, 0))
heisenberg <- within(heisenberg, site_dumb07<-ifelse(SITE=="PITT", 1, 0))
heisenberg <- within(heisenberg, site_dumb08<-ifelse(SITE=="SDSU", 1, 0))
heisenberg <- within(heisenberg, site_dumb09<-ifelse(SITE=="TRINITY", 1, 0))
heisenberg <- within(heisenberg, site_dumb10<-ifelse(SITE=="UM_1", 1, 0))
heisenberg <- within(heisenberg, site_dumb11<-ifelse(SITE=="UM_2", 1, 0))
heisenberg <- within(heisenberg, site_dumb12<-ifelse(SITE=="USM", 1, 0))
heisenberg <- within(heisenberg, site_dumb13<-ifelse(SITE=="YALE", 1, 0))
# 
# heisenberg[,"site_dumb01"] + heisenberg[,"site_dumb02"] + heisenberg[,"site_dumb03"] + heisenberg[,"site_dumb04"] + heisenberg[,"site_dumb05"] + heisenberg[,"site_dumb06"] + heisenberg[,"site_dumb07"] + heisenberg[,"site_dumb08"] + heisenberg[,"site_dumb09"] + heisenberg[,"site_dumb10"] + heisenberg[,"site_dumb11"] + heisenberg[,"site_dumb12"] +heisenberg[,"site_dumb13"]

#Correlation <- matrix(scan("/Volumes/Purkinje/SVM-RFE/ABIDE_allsubjects.txt", n = 252*24091 ), 252, 24091, byrow = TRUE)
# Correlation = allsubj
number_columns = ncol(regressed)
regressed.site2 <- mat.or.vec(nsubj, number_columns)
# regressed[,1] = allsubj$response
# 1 =ASD= <2>; 0 =TD= <1>
regressed.age = regressed
 
for( i in 1: number_columns) {  # start with 2nd column 

lm.r = lm( regressed.age[,i] ~ site_dumb01 + site_dumb02 + site_dumb03 + site_dumb04 + site_dumb05 + site_dumb06 + site_dumb07 + site_dumb08 + site_dumb09 + site_dumb10 + site_dumb11 + site_dumb12 + site_dumb13-1, data = heisenberg)
row <- resid(lm.r)	
 col <- t(row) 
 regressed.site2[,i] <- t(col)
print(i) 
}

model.matrix(lm.r)
summary(lm.r)



#
par(mfrow=c(1,1))

avg.abide = NULL
avg.age=NULL
avg.site=NULL
for (i in 1:nsubj){
  # avg.abide[i] = mean(as.matrix(abide449[i,]))
  avg.age[i]= mean(regressed.age[i,])
  avg.site[i]=mean(regressed.site2[i,])
}

boxplot(avg.abide ~ demo$SITE, data=demo)
boxplot(avg.age ~ demo$SITE, data=demo)
boxplot(avg.site ~ demo$SITE, data=demo)




write.csv(regressed, file = "regressed_data.csv")
