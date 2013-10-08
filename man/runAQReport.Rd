\name{runAQReport}
\alias{runAQReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{runAQReport - runing the AmpliQueso report
%%  ~~function to do ... ~~
}
\description{
This function creates AmpliQueso report in PDF format using either Article or Beamer LaTeX document classes.
Report is created in current working directory. Please make sure that you have Beamer LaTeX
package installed prior to running this report if you choose this output report type.
}
\usage{
runAQReport(iCovdesc, iBedFile, iRefSeqFile, iGroup, iT1, iT2, iTopN = 5, iMinQual, iReportFormat = "pdf", iReportType = "article",iReportPath="", iVerbose = TRUE, iVCF = NULL, iParallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iCovdesc}{
Covdesc-like data frame - BAM files are read from row names. Similar to covdesc in simpleaffy and rnaSeqMap - tab delimited table of BAM files and groupings. This is the table of experimental design.
}
  \item{iBedFile}{
BED file with same genomic coordinates as in BAMs (preferably standard hg19)
}
  \item{iRefSeqFile}{
Reference sequence in FASTA format used for SNPs detection. If set to NULL (by default) will skip SNP report.
}
  \item{iGroup}{
Name of the attribute/group in the experimental description (i.e. covdesc)
}
  \item{iT1}{
name of the first group from the experiment (it should be consistent with the groupings specified in covdesc file)
}
  \item{iT2}{
name of the other group from the experiment (it should be consistent with the groupings specified in covdesc file)
}
  \item{iTopN}{
number specifing top n expressed genes to be reported
}
  \item{iMinQual}{
number specifing the minimum quality of the SNPs returned by samtools used to filter out these of poor quality. Please refer to
      Samtools documentation for further details on the definition of this measure.
}
  \item{iReportFormat}{
Report output format. Currently PDF only.
}
  \item{iReportType}{
LaTeX document class to be used for report generation.
}
  \item{iReportPath}{
Output path for the report. It is created in the current working directory by default.
}
  \item{iVerbose}{
Boolean specifing whether R code should be inlcude in the report, if set to FALSE only figures and tables will be included.
}
  \item{iVCF}{
TBD
}
 \item{iParallel}{
Boolean specifing whether R code should be run in parallel -- it is turned on by default..
}

}

\value{
No vaue returned.
}

\author{
   Alicja Szabelska, Marek Wiewiorka, Michal Okoniewski
}
\examples{
#runAQReport(iCovdesc="covdesc",iBedFile="test.BED",iRefSeqFile=NULL,iGroup="group",iT1="sick",iT2="healthy",iTopN=5,iMinQual=10,iReportFormat="pdf",iReportType="article",iVerbose=FALSE,iVCF=NULL)
}
