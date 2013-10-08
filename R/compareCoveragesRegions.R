#includes
compareCoveragesReg <- function(iBedFile, iGroup, iT1, iT2,  iMeasure=c("DA","QQ", "PP", "HD1", "HD2"), iCovdesc="covdesc",iParallel=TRUE){
  if(iParallel == TRUE)
     #assume that the number of regions is greatet than measures which is true in most cases
    .compareCoveragesRegParR(iBedFile, iGroup, iT1, iT2,iCovdesc)
  else
    .compareCoveragesRegSeq(iBedFile, iGroup, iT1, iT2,  iMeasure,iCovdesc)
  
}



.compareCoveragesRegSeq <- function(iBedFile, iGroup, iT1, iT2,  iMeasure=c("DA","QQ", "PP", "HD1", "HD2"), iCovdesc="covdesc"){
	#read BED file
	regions <- read.table(iBedFile,col.names=c("chr", "st", "en","V4","V5","name")) 
	covResult <- matrix(nrow=nrow(regions), ncol=length(iMeasure) )
	#colnames(covResult) <- iMeasures
	#rownames(covResult) <- as.character(regions$name)
	#iterate over measures	
	for(m in 1:length(iMeasure) ){
     
	   #interate over regions
	   for ( i in 1:nrow(regions) ){
		
	     	region <- regions[i,]
	     	chrName <- as.character(region$chr)
	     	startPos <- region$st
	     	endPos   <- region$en
	    	covResult[i,m] <-compareCoverages(ch=chrName, st=startPos, en=endPos, group=iGroup, iT1, iT2, measure=iMeasure[m], covdesc=iCovdesc)	
	  }
	
	}
	return (covResult)	
}

##### profiling and performance testing

.compareCoveragesRegApply <- function(iBedFile, iGroup, iT1, iT2,  iMeasure=c("DA","QQ", "PP", "HD1", "HD2"), iCovdesc="covdesc"){
        #read BED file
        regions <- read.table(iBedFile,col.names=c("chr", "st", "en","V4","V5","name"))
        covResult <- matrix(nrow=nrow(regions), ncol=length(iMeasure) )
        #colnames(covResult) <- iMeasures
        #rownames(covResult) <- as.character(regions$name)
        #iterate over measures
        m=NULL
        for(m in 1:length(iMeasure) ){

          
          covResult[,m] <- apply(regions,1, .callCompareCoverages, iGroup=iGroup, iT1=iT1, iT2=iT2, iMeasure=iMeasure[m],iCovdesc=iCovdesc) 
	  #interate over regions
           #for ( i in 1:nrow(regions) ){

          #      region <- regions[i,]
          #      chrName <- as.character(region$chr)
          #      startPos <- region$st
          #      endPos   <- region$en
          #      covResult[i,m] <-compareCoverages(ch=chrName, st=startPos, en=endPos, group=iGroup, iT1, iT2, measure=iMeasure[m], covdesc=iCovdesc)
        	 
	# }

        }
        return (covResult)
}



.callCompareCoverages <- function (region, iGroup, iT1, iT2, iMeasure, iCovdesc){

  chrName <- as.character(region[1])
	startPos <- region[2]
  endPos   <- region[3]
  result <-compareCoverages(ch=chrName, st=startPos, en=endPos, group=iGroup, iT1, iT2, measure=iMeasure, covdesc=iCovdesc)
	return (result)
}

.callCompareCoveragesBatch <- function (region, iGroup, iT1, iT2, iCovdesc){
  
  chrName <- as.character(region[1])
  startPos <- region[2]
  endPos   <- region[3]
  result <-compareCoverages(ch=chrName, st=startPos, en=endPos, group=iGroup, iT1, iT2,measure=c("DA","QQ", "PP", "HD1", "HD2"), covdesc=iCovdesc)
  return (result)
}


#parallel by measure
.compareCoveragesRegParM <- function(iBedFile, iGroup, iT1, iT2,  iMeasure=c("DA","QQ", "PP", "HD1", "HD2"), iCovdesc="covdesc"){

  .loadLibs("parallel")
  .loadLibs("doParallel")
  .loadLibs("foreach")

  regions <- read.table(iBedFile,col.names=c("chr", "st", "en","V4","V5","name"))
  cl <- makeCluster(detectCores())
  #manual detection
  #cl <- makeCluster(6)
  registerDoParallel(cl)
  covResult <- matrix(nrow=nrow(regions), ncol=length(iMeasure) )
  m=NULL
  covResult <-foreach( m=1:length(iMeasure), .export=c(".compareCoveragesRegApply","iMeasure","iCovdesc", "iT1", "iT2", "iGroup", "iBedFile"),.combine=cbind, .packages="ampliQueso" ) %dopar% {
     #icallSamPileup(i, covdesc,minQual,refSeqFile, bedFile  )
     .compareCoveragesRegApply(iBedFile, iGroup, iT1, iT2, iMeasure[m], iCovdesc)
  	
	} 

  stopCluster(cl)	
  return (covResult)
}


#parallel by region uses compareCoveragesBatch instead of compareCoverages
.compareCoveragesRegParR <- function(iBedFile, iGroup, iT1, iT2, iCovdesc="covdesc"){
  
  .loadLibs("parallel")
  .loadLibs("doParallel")
  .loadLibs("foreach")
  regions <- read.table(iBedFile,col.names=c("chr", "st", "en","V4","V5","name"))
  cl <- makeCluster(detectCores())
  #manual detection
  #cl <- makeCluster(6)
  registerDoParallel(cl)
  #temp hardcode for ncol = number of measures
  covResult <- matrix(nrow=nrow(regions), ncol=5)
  r=NULL
  covResult <-foreach( r=1:nrow(regions), .export=c(".callCompareCoveragesBatch","compareCoverages", "iBedFile"),.combine=rbind, .packages="ampliQueso" ) %dopar% {
    .callCompareCoveragesBatch(as.matrix(regions[r,]), iGroup, iT1, iT2, iCovdesc)
    
  } 
  
  stopCluster(cl)	
  return (covResult)
}