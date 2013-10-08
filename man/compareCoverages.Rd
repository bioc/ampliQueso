\name{compareCoverages}
\alias{compareCoverages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{compareCoverages - camel measure for a single comparison.
}

\description{
Comparing coverage profiles with the camel measures for the two groups of BAM files, for a given genomic region. 
}

\usage{
compareCoverages(ch, st, en, group, t1, t2, measure="DA", covdesc="covdesc")
}
\arguments{
      \item{ch}{Chromosome of the region of interest 
    }
      \item{st}{Start of the region of interest
    }
      \item{en}{End of the region of interest
    }
      \item{group}{Name of the attribute/group in the experimental description (i.e. covdesc)
    }
     \item{t1}{First name of the attribute used for grouping
    }
     \item{t2}{Another name of the attribute used for grouping
    }
      \item{measure}{Method of comparison for the coverages in the selected genomic region
    }
      \item{covdesc}{Covdesc-like data frame - BAM files are read from row names
    }
}

\value{
A single camel value with the results of comparison.
}
\author{
   Alicja Szabelska, Marek Wiewiorka, Michal Okoniewski
}



