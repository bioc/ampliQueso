# compares two NuclDistributions according to camel measures
compareCoverages <- function(ch, st, en, group, t1, t2, measure="DA", covdesc="covdesc")
{
  # strand hardcoded, as we assume non-stranded protocol
  str <- 1
  cvd<-read.table(covdesc)
  numbams <- dim(cvd)[1]
  samples1 <- which(cvd[,group]==t1)
  samples2 <- which(cvd[,group]==t2)
  rs  <- newSeqReads(ch,st, en, str)
  rs <- getBamData(rs,1:numbams, covdesc.file=covdesc)
  nd <- getCoverageFromRS(rs,1:numbams)
  v1 <- averageND(nd, samples1)
  v2 <- averageND(nd, samples2)
  v1Data <- getData(v1)
  v2Data <- getData(v2)
  dd <- cbind(RleList2matrix(v1Data), RleList2matrix(v2Data))
  #dd <- cbind(RleList2matrix(v1@data), RleList2matrix(v2@data))
  out <- c()
  if ("DA" %in% measure ) out <- cbind(out,diff_area(dd))
  if ("QQ" %in% measure ) out <- cbind(out,qq_plot(dd))
  if ("PP" %in% measure ) out <- cbind(out,pp_plot(dd))
  if ("HD1" %in% measure ) out <- cbind(out,hump_diff1(dd))
  if ("HD2" %in% measure ) out <- cbind(out,hump_diff2(dd))
  out
}

# compareCoveragesBatch <- function(ch, st, en, group, t1, t2, covdesc="covdesc")
# {
#   # strand hardcoded, as we assume non-stranded protocol
#   str <- 1
#   cvd<-read.table(covdesc)
#   numbams <- dim(cvd)[1]
#   samples1 <- which(cvd[,group]==t1)
#   samples2 <- which(cvd[,group]==t2)
#   rs  <- newSeqReads(ch,st, en, str)
#   rs <- getBamData(rs,1:numbams, covdesc=covdesc)
#   nd <- getCoverageFromRS(rs,1:numbams)
#   v1 <- averageND(nd, samples1)
#   v2 <- averageND(nd, samples2)
#   dd <- cbind(RleList2matrix(v1@data), RleList2matrix(v2@data))
#   out <- diff_area(dd)
#   out <- cbind(out,qq_plot(dd))
#   out <- cbind(out,pp_plot(dd))
#   out <- cbind(out,hump_diff1(dd))
#   out <- cbind(out,hump_diff2(dd))
#   out
# }

