design <-read.table("DATA_allpoints",header=FALSE, sep = " ")
design <- as.matrix(design[c(12, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 42, 43, 44, 45, 47, 48, 51, 52, 53, 54, 55, 60, 61, 63, 64, 65, 66, 67, 68, 69, 70, 72, 73, 76, 78, 79, 81, 82, 83, 85),])
rownames(design) <- c()
write.table(design,file = "DATA_filtered",sep = " ",col.names = FALSE,row.names = FALSE)
