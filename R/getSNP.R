getSNP <-function(covdesc,minQual,refSeqFile, bedFile=NULL,iParallel=TRUE){
  
  if(iParallel == TRUE)
    .getSNPPar(covdesc,minQual,refSeqFile, bedFile)
  else
    .getSNP(covdesc,minQual,refSeqFile, bedFile)
  
  
}


.callSamPileup <-function(procId, covdesc,minQual,refSeqFile, bedFile){
	
  #regions limited by BED file
  if( hasArg(bedFile) && !is.null(bedFile)){
    system(paste("samtools mpileup -uf ",refSeqFile,rownames(covdesc)[procId],"-l", bedFile, 
                 paste("| bcftools view -vcg - > /tmp/aq_tmp_", procId, ".vcf",sep="") , sep=" " ) )
  }  
  else{
    system(paste("samtools mpileup -uf ",refSeqFile,rownames(covdesc)[procId], "-l null", 
                 paste("| bcftools view -vcg - > /tmp/aq_tmp_", procId, ".vcf",sep="") , sep=" " ) )
    
    
  }
    aqVcf <- readVcf(paste("/tmp/aq_tmp_", procId, ".vcf", sep=""), "hg19")
    aqDF<-data.frame(sequence=as.character(seqnames(aqVcf)),start=start(ranges(aqVcf)),
                     end=end(ranges(aqVcf)),width=width(ranges(aqVcf)),ref=as.character(ref(aqVcf)),
                     strand=as.character(strand(aqVcf)),quality=qual(aqVcf) )
    aqSNPDF <- aqDF[aqDF$quality > minQual,]
    return (aqSNPDF)

}


.getSNP <- function(covdesc,minQual,refSeqFile, bedFile=NULL){
#load libs
  .loadLibs("VariantAnnotation")
#init global data structures
###
  aqSNPL=list()
  #load covdesc file  
  covdescT<-read.table(covdesc,header=TRUE)
  for ( i in 1:nrow(covdescT) ) {
    #regions limited by BED file
    if( hasArg(bedFile) && !is.null(bedFile)){
        system(paste("samtools mpileup -uf ",refSeqFile,rownames(covdescT)[i],"-l", bedFile, 
                    "| bcftools view -vcg - > /tmp/aq_tmp.vcf" , sep=" " ) )
    }
    else{
      system(paste("samtools mpileup -uf ",refSeqFile,rownames(covdescT)[i], "-l null",
                   "| bcftools view -vcg - > /tmp/aq_tmp.vcf" , sep=" " ) )
      
    }
    
    aqVcf <- readVcf("/tmp/aq_tmp.vcf", "hg19")
    aqDF<-data.frame(sequence=as.character(seqnames(aqVcf)),start=start(ranges(aqVcf)),
                     end=end(ranges(aqVcf)),width=width(ranges(aqVcf)),ref=as.character(ref(aqVcf)),
                     strand=as.character(strand(aqVcf)),quality=qual(aqVcf) )
    aqSNPL[[length(aqSNPL)+1]] <- aqDF[aqDF$quality > minQual,]
    #covdescT$group[i]
  #print(test[i] ) 
    
  }
  return (aqSNPL)

}
.getSNPPar <- function(covdesc,minQual,refSeqFile, bedFile=NULL){
  .loadLibs("parallel")
  .loadLibs("doParallel")
  .loadLibs("foreach")
  .loadLibs("VariantAnnotation")
  #automatic no of cores detection
  cl <- makeCluster(detectCores())
  #manual detection
  #cl <- makeCluster(6)
  registerDoParallel(cl)
  aqSNPL=list()
  covdescT<-read.table(covdesc,header=TRUE)
  legendPaths<-strsplit( rownames(covdescT),"/")
  i=NULL
  aqSNPL<- foreach( i=1:nrow(covdescT), .export=c(".callSamPileup","covdescT"), .packages=c("VariantAnnotation","ampliQueso")) %dopar% {
     .callSamPileup(i, covdescT,minQual,refSeqFile, bedFile  )
  }
  
  stopCluster(cl)
 
  
  return (aqSNPL)
}

