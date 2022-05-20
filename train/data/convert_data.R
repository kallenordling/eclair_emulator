design <-read.csv("design.csv",header=TRUE)
design <- design[,-1]

zeros <- rep(0,nrow(design))

response <-read.csv("cfrac.csv",header=FALSE)

out <- cbind(design,zeros,response)

write.table(out,file = "DATA",sep = " ",col.names = FALSE,row.names = FALSE)