\name{ndMin}
\alias{ndMin}
\docType{data}
\title{
ampliQueso sample data
}
\description{
Nucleotide distribution object for a region  with the lowest camel measures difference between samples.
}
\usage{data(ampliQueso)}
\format{
  The format is:
Formal class 'NucleotideDistr' [package "rnaSeqMap"] with 13 slots
  ..@ chr              : chr "chr2"
  ..@ start            : num 68615536
  ..@ end              : num 68620354
  ..@ strand           : num 1
  ..@ type             : chr(0) 
  ..@ data             :List of 8
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:78] 0 1 47 48 49 54 500 503 504 505 ...
  .. .. .. ..@ lengths        : int [1:78] 12 1 1 9 1 1 1 3 2 6 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:87] 0 2 61 786 787 788 790 793 795 796 ...
  .. .. .. ..@ lengths        : int [1:87] 12 1 12 1 1 1 1 1 2 3 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:72] 0 1 48 49 50 52 506 508 511 513 ...
  .. .. .. ..@ lengths        : int [1:72] 12 1 9 1 1 1 2 2 1 3 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:74] 0 60 62 63 65 738 741 743 746 748 ...
  .. .. .. ..@ lengths        : int [1:74] 13 1 7 3 1 1 2 2 4 2 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:61] 0 20 23 513 516 517 518 519 520 521 ...
  .. .. .. ..@ lengths        : int [1:61] 13 11 1 1 3 1 1 1 3 2 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:61] 0 19 21 437 438 440 441 442 444 445 ...
  .. .. .. ..@ lengths        : int [1:61] 13 11 1 1 2 1 1 3 7 5 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:64] 0 1 36 37 466 467 468 469 470 472 ...
  .. .. .. ..@ lengths        : int [1:64] 12 1 2 10 1 1 1 1 2 1 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:84] 0 1 51 52 53 54 725 726 728 729 ...
  .. .. .. ..@ lengths        : int [1:84] 12 1 5 3 3 1 1 3 3 1 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  ..@ assayData        :<environment: 0x104fc118> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	0 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	0 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  .. .. ..@ name             : chr ""
  .. .. ..@ lab              : chr ""
  .. .. ..@ contact          : chr ""
  .. .. ..@ title            : chr ""
  .. .. ..@ abstract         : chr ""
  .. .. ..@ url              : chr ""
  .. .. ..@ pubMedIds        : chr ""
  .. .. ..@ samples          : list()
  .. .. ..@ hybridizations   : list()
  .. .. ..@ normControls     : list()
  .. .. ..@ preprocessing    : list()
  .. .. ..@ other            : list()
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ annotation       : chr(0) 
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	0 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 3
  .. .. .. ..$ : int [1:3] 3 0 1
  .. .. .. ..$ : int [1:3] 2 20 1
  .. .. .. ..$ : int [1:3] 1 3 0
}
\examples{
data(ampliQueso)
}
\keyword{datasets}
