\name{compareCoveragesReg}
\alias{compareCoveragesReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{compareCoveragesReg - calculating given camel measures for a list of regions specified by a BED file (sequential and parallel versions)
%%  ~~function to do ... ~~
}
\description{
 Calculating given camel measures for a list of regions specified by a BED file (sequential version)}
\usage{
compareCoveragesReg(iBedFile, iGroup, iT1, iT2, iMeasure = c("DA", "QQ", "PP", "HD1", "HD2"), iCovdesc = "covdesc",iParallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
      \item{iBedFile}{
     BED file with same genomic coordinates as in BAMs (preferably standard hg19)
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
      \item{iMeasure}{
    vector specifing camel measures to be calculated
    }
      \item{iCovdesc}{
    Covdesc-like data frame - BAM files are read from row names. Similar to covdesc in simpleaffy and rnaSeqMap - tab delimited table of BAM files and groupings. This is the table of experimental design.
    }
     \item{iParallel}{
    boolean specifing whether calculations shoulb be run in parallel,true by default,
    The degree of parallelism is set up automagically and by default equals to the number of logical CPU cores.
    }
}

\value{
data frame with camel measures calculated for the specified regions and samples}
\author{
   Alicja Szabelska, Marek Wiewiorka, Michal Okoniewski
}


