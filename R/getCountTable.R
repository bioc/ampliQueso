getCountTable <- function(covdesc="covdesc", bedFile="amplicons.bed")
{
    tout <- NULL
    annot<-read.table(bedFile)
    numgenes <- dim(annot)[1]
    cvd <- read.table(covdesc)
    numbams <- dim(cvd)[1]
    for (i in 1:numgenes)
    {
        ch <- as.character(annot[i,1]);
        st <- as.numeric(annot[i,2]);
        en <- as.numeric(annot[i,3]);
        str <- as.character(annot[i,4]);
        str <- .strand2number(str)
        rs <- rnaSeqMap:::newSeqReads(ch,st, en, str);
        rs <- getBamData(rs,1:numbams, covdesc.file=covdesc)
        rsData <- getData(rs)
	      counts<-sapply(rsData, length)
      	#if (normalize=="length") counts <- normLen(unlist(counts), rs, c=normLenConst)
	      tout <- rbind(tout, counts)
        cat('.')
    }
    rownames(tout) <- annot[,6]
    colnames(tout) <- rownames(cvd)
    tout
}


.strand2number <- function(s)
{
  out <- NULL
  if (identical(s,"1") || identical(s,"+") ) out <- 1
  else if (identical( s,"-1") || identical(s,"-") ) out <- (-1)
  else out <- 1
  out
}
