\name{camelTest}
\alias{camelTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{camelTest - calculating permutation tests for camel measures (sequential and parallel versions)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculating permutation tests for reads coverage (camel) measures for regions specified by a BED file.
This function features a few normalization modes and camel measures.
}
\usage{
camelTest(iBedFile, iCovdesc = "covdesc", iT1, iT2, iNorm = c("none", "density", "minMax"), iMeasure = c("DA", "QQ", "PP", "HD1", "HD2"), iSizes = NULL, iNPerm = 100,iParallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
      \item{iBedFile}{
    BED file with genomic coordinates refering to the same genome as the BAMs
    }
      \item{iCovdesc}{Covdesc-like file  - BAM files are read from row names. Similar to covdesc in simpleaffy and rnaSeqMap - tab delimited table of BAM files and groupings. This is the table of experimental design. 
    }
      \item{iT1}{
    name of the first group from the experiment (it should be consistent with the groupings specified in covdesc file)
    }
      \item{iT2}{
    name of the other group from the experiment (it should be consistent with the groupings specified in covdesc file)
    }
      \item{iNorm}{
    vector specifing the data normalisations prior to camel measures calculations
    }
      \item{iMeasure}{
     camel measure to be used in the non-parametric test. Camel measures are implemented in the package rnaSeqMap and described in Okoniewski et al, NAR, 2011}
      \item{iSizes}{
    TBD}
    \item{iNPerm}{
    integer number specifing the number of permutations for camel tests calcualtions}
    
    \item{iParallel}{
    boolean specifing whether calculations shoulb be run in parallel,true by default,
    The degree of parallelism is set up automagically and by default equals to the number of logical CPU cores.
    }
}
\value{
data frame with p-values of camel tests for all regions, measures and data normalisations specified}
\author{
    Alicja Szabelska, Marek Wiewiorka, Michal Okoniewski
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
