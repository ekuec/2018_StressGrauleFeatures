library(pROC)
library(MASS)
library(PRROC)

data <- read.table("human_pme.data")
sg   <- read.table("nsg_human.txt")
cty  <- read.table("human_ncyto.txt")
nuc  <- read.table("human_nucleus.txt")
nol  <- read.table("human_nucleolus.txt")
rnp  <- read.table("human_RNPgranule.txt")
pby  <- read.table("human_pbody.txt")
llps <- read.table("human_llps.txt")
ggs  <- read.table("gingras_gold.txt")

names(data) <- c("name","diso","abd","csl","int","len","max","phs","pip","rna","mrf","lps","cat","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
names(sg)  <- c("name")
names(ggs) <- c("name")

data$Sabd=(log(data$abd,10)-mean(log(data$abd,10)))/sd(log(data$abd,10))
data$Scsl=(data$csl-mean(data$csl))/sd(data$csl)
data$Sdiso=(data$diso-mean(data$diso))/sd(data$diso)
data$Slen=(log(data$len,10)-mean(log(data$len,10)))/sd(log(data$len,10))
data$Spip=(data$pip-mean(data$pip))/sd(data$pip)
data$Sint=(data$int-mean(data$int))/sd(data$int)
data$Sphs=(data$phs-mean(data$phs))/sd(data$phs)
data$Smrf=(data$mrf-mean(data$mrf))/sd(data$mrf)
data$Slps=(data$lps-mean(data$lps))/sd(data$lps)
data$Smax=(data$max-mean(data$max))/sd(data$max)
data$SA=(data$A-mean(data$A))/sd(data$A)
data$SC=(data$C-mean(data$C))/sd(data$C)
data$SD=(data$D-mean(data$D))/sd(data$D)
data$SE=(data$E-mean(data$E))/sd(data$E)
data$SF=(data$F-mean(data$F))/sd(data$F)
data$SG=(data$G-mean(data$G))/sd(data$G)
data$SH=(data$H-mean(data$H))/sd(data$H)
data$SI=(data$I-mean(data$I))/sd(data$I)
data$SK=(data$K-mean(data$K))/sd(data$K)
data$SL=(data$L-mean(data$L))/sd(data$L)
data$SM=(data$M-mean(data$M))/sd(data$M)
data$SN=(data$N-mean(data$N))/sd(data$N)
data$SP=(data$P-mean(data$P))/sd(data$P)
data$SQ=(data$Q-mean(data$Q))/sd(data$Q)
data$SR=(data$R-mean(data$R))/sd(data$R)
data$SS=(data$S-mean(data$S))/sd(data$S)
data$ST=(data$T-mean(data$T))/sd(data$T)
data$SV=(data$V-mean(data$V))/sd(data$V)
data$SW=(data$W-mean(data$W))/sd(data$W)
data$SY=(data$Y-mean(data$Y))/sd(data$Y)

data$sg    <- as.numeric(data$name %in% sg$name)
data$ggs   <- as.numeric(data$name %in% ggs$name)
data$valid <- as.numeric(data$sg==0 & data$ggs==1)

folds = 10
repeats = 1

my.min.auc = c()
my.min.multi = c()
my.min.model = c()
my.min.zs = c()

for ( bigloop in 1:repeats ) {
  my.pool          <- data[which(data$valid != 1),]
  my.data          <- data[sample(nrow(my.pool)),]
  my.type          <- cut(seq(1,nrow(my.data)),breaks=2,labels=F)
  my.training      <- my.data[which(my.type==1,arr.ind=T),]
  my.control       <- my.data[which(my.type==2,arr.ind=T),]
  my.out           <- data[,c("name","sg")]
  my.min.out       <- data[,c("name","sg")]
  my.raw           <- data[,c("name","sg")]
  my.min.raw       <- data[,c("name","sg")]
  my.folds         <- cut(seq(1,nrow(my.training)),breaks=folds,labels=F)

 my.print <- my.training
 write.table(my.print, file='mags_set_test.tsv', quote=FALSE, sep='\t', col.names = NA)
 my.print <- my.control
 write.table(my.print, file='mags_set_train.tsv', quote=FALSE, sep='\t', col.names = NA)
 my.print <- data[which(data$valid == 1),]
 write.table(my.print, file='mags_set_valid.tsv', quote=FALSE, sep='\t', col.names = NA)

  
  for ( i in 1:folds) {
    my.test  <- my.training[which(my.folds==i,arr.ind=T),]
    my.train <- my.training[-which(my.folds==i,arr.ind=T),]
  
    my.min.mod <- glm(sg ~
           Sabd*(Sdiso +
           Scsl +
           Sphs +
           Spip +
           rna  +
           SL   +
           SG   )
           ,data=my.train,family="quasibinomial")
  
    my.min.roc <- roc(data$sg,predict(my.min.mod,data))
    my.min.auc <- rbind(my.min.auc, my.min.roc$auc)
    
    j=paste("mod",i,sep="")
    x=predict(my.min.mod,data)
    my.min.raw[,j] <- x
    my.min.out[,j] <- signif((x-mean(x))/sd(x),3)
  }

  my.min.out$AVG  <-rowMeans(my.min.out[,3:ncol(my.min.out)])
  my.min.raw$AVG  <-rowMeans(my.min.raw[,3:ncol(my.min.raw)])
  
  my.min.multi <- cbind(my.min.multi, my.min.raw$AVG)
  my.min.zs    <- cbind(my.min.zs,my.min.out$AVG)
  my.min.model <- cbind(my.min.model, my.min.mod$coefficients)

}

# ROC CURVES #
mycol=c("darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4",
        "darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4",
        "darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4",
        "darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4",
        "darkolivegreen","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4")
pdf("mags_roc.pdf")
plot(1,type='n',xlab="Specificity",ylab="Sensitivity",xlim=c(1,0),ylim=c(0,1),main="ROC Curves")
for ( i  in 1:ncol(my.min.multi) ) {
 lines(roc(data$sg,my.min.multi[,i]),lwd=2,col=mycol[i]) 
}
n1 <- paste("Model:",signif(mean(my.min.auc),4),sep=" ")
n1 <- paste(n1, "+/-", sep=" ")
n1 <- paste(n1, signif(sd(my.min.auc),2),sep=" ")
lines(roc(data$sg,data$cat),col="azure4",lty=1,lwd=3)
lines(roc(data$sg,data$pip),col="azure3",lty=1,lwd=3)
abline(1,-1,col="black",lwd=1,lty=2)
n2 <- paste("catGranule:",signif(roc(data$sg,data$cat)$auc,4),sep=" ")
n3 <- paste("PScore:",signif(roc(data$sg,data$pip)$auc,4),sep=" ")
legend("bottomright",legend=c(n1,n2,n3), col=c("darkolivegreen4","azure4","azure3"), lty=c(1,1,1),lwd=3)
dev.off()

# PRINT LISTS #
my.min.multi <- as.data.frame(my.min.multi)
my.min.zs    <- as.data.frame(my.min.zs)
my.min.model <- as.data.frame(my.min.model)

my.min.multi$AVG <- my.min.multi[,1]
my.min.zs$AVG    <- my.min.zs[,1]
my.min.model$AVG <- my.min.model[,1]

my.min.multi$name <- data$name
my.min.zs$name    <- data$name 

my.print <- my.min.multi[order(-my.min.multi$AVG),] 
write.table(my.print, file='mags_rawscore.tsv', quote=FALSE, sep='\t', col.names = NA)
my.print <- my.min.zs[order(-my.min.zs$AVG),]
write.table(my.print, file='mags_zscore.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(my.min.model, file='out_model.tsv', quote=FALSE, sep='\t', col.names = NA)

# PR CURVES #
my.min.multi$sg<-data$sg
my.min.multi$cat<-data$cat
my.min.multi$pip<-data$pip
pdf("mags_pr.pdf")
plot(1,type='n',xlab="Fraction Recall",ylab="Precision",xlim=c(0,1),ylim=c(0,1),main="")
tmp.pr <- pr.curve(scores.class0 = my.min.multi[which(my.min.multi$sg==1),"AVG"], scores.class1 = my.min.multi[which(my.min.multi$sg==0),"AVG"],curve=T)
plot(tmp.pr,col="darkolivegreen4",lwd=3,add=T)
tmp.pr <- pr.curve(scores.class0 = my.min.multi[which(my.min.multi$sg==1),"pip"], scores.class1 = my.min.multi[which(my.min.multi$sg==0),"pip"],curve=T)
plot(tmp.pr,col="azure3",lwd=3,add=T)
tmp.pr <- pr.curve(scores.class0 = my.min.multi[which(my.min.multi$sg==1),"cat"], scores.class1 = my.min.multi[which(my.min.multi$sg==0),"cat"],curve=T)
plot(tmp.pr,col="azure4",lwd=3,add=T)
legend("topright",legend=c("Model","catGranule","PScore"), col=c("darkolivegreen4","azure4","azure3"),lty=1,lwd=3)
dev.off()
