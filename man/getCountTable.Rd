\name{getCountTable}
\alias{getCountTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getCountTable - obtaining typical count table.
}

\description{
Obtaining a count table for all the amplicons described in the BED file. 
}

\usage{
getCountTable(covdesc="covdesc", bedFile="amplicons.bed")
}
\arguments{
  \item{covdesc}{Covdesc-like data frame - BAM files are read from row names. Similar to covdesc in simpleaffy and rnaSeqMap - tab delimited table of BAM files and groupings. This is the table of experimental design. 
}
  \item{bedFile}{BED file with same genomic coordinates as in BAMs (preferably standard hg19 )
}
}

\value{
Count table with the number of columns the same as number of samples and the number of rows the same as number of amplicons.
}
\author{
   Alicja Szabelska, Marek Wiewiorka, Michal Okoniewski
}
\examples{
library(ampliQueso)
setwd(path.package("ampliQueso"))
cc <- getCountTable(covdesc=system.file("extdata","covdesc",package="ampliQueso") , bedFile=system.file("extdata","AQ.bed",package="ampliQueso"))

}
