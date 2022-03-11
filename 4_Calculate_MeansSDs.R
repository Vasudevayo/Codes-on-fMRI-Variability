#Written by Jeffrey Malins

for(mask in c("LIFGop","LIFGpt","LITG","LMTG","LSPL","LSTG","LThal")){

#Load in data from ROIs
#Speech
MeanA <- read.table(paste(mask,"_Means_AudWordMM_Correct.txt",sep=""), header = FALSE, fill=TRUE, col.names=c("Subj",paste0('Trial', seq_len(40))))
MeanA_noSubj <- as.matrix(MeanA[,-1])
OutliersA <- read.table(paste(mask,"_Outliers_AudWordMM_Correct.txt",sep=""), header = FALSE, fill=TRUE, col.names=paste0('Trial', seq_len(40)))
OutliersA_noSubj <- as.matrix(OutliersA[,-1])
mean_outliersA <- apply(OutliersA_noSubj,1,mean,na.rm=TRUE)
meansA <- apply(MeanA_noSubj,1,mean,na.rm=TRUE)
sdsA <- apply(MeanA_noSubj,1,bootstrap::jackknife,sd,na.rm=TRUE)
sd_jackA <- array(,dim=nrow(MeanA_noSubj))
for(i in 1:nrow(MeanA_noSubj)){sd_jackA[i]<-mean(sdsA[[i]]$jack.values)}

#Print
MeanV <- read.table(paste(mask,"_Means_VisWordMM_Correct.txt",sep=""), header = FALSE, fill=TRUE, col.names=paste0('Trial', seq_len(40)))
MeanV_noSubj <- as.matrix(MeanV[,-1])
OutliersV <- read.table(paste(mask,"_Outliers_VisWordMM_Correct.txt",sep=""), header = FALSE, fill=TRUE, col.names=paste0('Trial', seq_len(40)))
OutliersV_noSubj <- as.matrix(OutliersV[,-1])
mean_outliersV <- apply(OutliersV_noSubj,1,mean,na.rm=TRUE)
meansV <- apply(MeanV_noSubj,1,mean,na.rm=TRUE)
sdsV <- apply(MeanV_noSubj,1,bootstrap::jackknife,sd,na.rm=TRUE)
sd_jackV <- array(,dim=nrow(MeanV_noSubj))
for(i in 1:nrow(MeanV_noSubj)){sd_jackV[i]<-mean(sdsV[[i]]$jack.values)}

#Write out files
Stats <- data.frame(MeanA[,1],meansA,sd_jackA,mean_outliersA,meansV,sd_jackV,mean_outliersV)
names(Stats) <- c("Subject_ID",paste(mask,"_meanA",sep=""),paste(mask,"_sdA",sep=""),paste(mask,"_OutliersA",sep=""),paste(mask,"_meanV",sep=""),paste(mask,"_sdV",sep=""),paste(mask,"_OutliersV",sep=""))
write.csv(Stats,paste(mask,"Stats.csv",sep=""), row.names = FALSE, quote = FALSE)

}

