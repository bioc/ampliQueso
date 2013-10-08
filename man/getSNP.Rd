\name{getSNP}
\alias{getSNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getSNP - finding SNPs in all samples for regions specified by a BED file all regions if BED file not provided (sequential and parallel versions)
%%  ~~function to do ... ~~
}
\description{
Finding SNPs in all samples for regions specified by a BED file (sequential version) using Samtools mpileup functionality.
Please make sure that samtools executables are installed and available from command line without specifing the full path(e.g. added to the PATH environment variable)
}
\usage{
getSNP(covdesc, minQual, refSeqFile, bedFile = NULL,iParallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
      \item{covdesc}{
    Covdesc-like data frame - BAM files are read from row names. Similar to covdesc in simpleaffy and rnaSeqMap - tab delimited table of BAM files and groupings. This is the table of experimental design.
    }
      \item{minQual}{
      number specifing the minimum quality of the SNPs returned by samtools used to filter out these of poor quality. Please refer to
      Samtools documentation for further details on the definition of this measure.
    }
      \item{refSeqFile}{
      Reference sequence in FASTA format used for SNPs detection.
    }
      \item{bedFile}{
    BED file with same genomic coordinates as in BAMs (preferably standard hg19)
    }
    \item{iParallel}{
    boolean specifing whether calculations shoulb be run in parallel,true by default,
    The degree of parallelism is set up automagically and by default equals to the number of logical CPU cores.
    }
}

\value{
list of data frames, for one sample each containing coordinates of the SNP with its quality
}

\author{
   Alicja Szabelska, Marek Wiewiorka, Michal Okoniewski
}
