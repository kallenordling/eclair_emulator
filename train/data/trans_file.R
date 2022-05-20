design <-read.csv("design.csv",header=TRUE)
design <- design[,-1]

meanvec <- colMeans(design)
stdvec <- apply(design,2,sd)

outvec <- rbind(meanvec,stdvec)

write.csv(outvec,file = "trans_data.csv")