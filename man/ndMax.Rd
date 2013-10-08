\name{ndMax}
\alias{ndMax}
\docType{data}
\title{
ampliQueso sample data
}
\description{
Nucleotide distribution object for a region  with the largest camel measures difference between samples.
}
\usage{data(ampliQueso)}
\format{
  The format is:
Formal class 'NucleotideDistr' [package "rnaSeqMap"] with 13 slots
  ..@ chr              : chr "chr1"
  ..@ start            : num 1.93e+08
  ..@ end              : num 1.93e+08
  ..@ strand           : num 1
  ..@ type             : chr(0) 
  ..@ data             :List of 8
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:68] 0 15 16 159 164 165 167 170 175 176 ...
  .. .. .. ..@ lengths        : int [1:68] 11 10 1 1 2 2 1 3 1 2 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:45] 0 4 5 6 77 80 81 82 83 84 ...
  .. .. .. ..@ lengths        : int [1:45] 11 7 3 1 1 2 1 8 3 1 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:78] 0 14 15 230 234 238 239 240 243 245 ...
  .. .. .. ..@ lengths        : int [1:78] 11 10 1 1 2 1 1 1 3 2 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:49] 0 1 11 12 95 97 98 100 101 103 ...
  .. .. .. ..@ lengths        : int [1:49] 10 1 10 1 3 1 2 7 1 2 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:64] 0 8 9 10 118 121 124 125 127 131 ...
  .. .. .. ..@ lengths        : int [1:64] 11 5 4 2 1 2 1 1 1 1 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:37] 0 4 47 48 49 50 51 52 53 55 ...
  .. .. .. ..@ lengths        : int [1:37] 11 11 1 3 10 2 6 1 1 1 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:77] 0 1 26 28 252 259 262 264 269 271 ...
  .. .. .. ..@ lengths        : int [1:77] 10 1 10 1 1 2 1 2 2 1 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'Rle' [package "IRanges"] with 4 slots
  .. .. .. ..@ values         : int [1:71] 0 16 167 171 174 175 177 179 180 181 ...
  .. .. .. ..@ lengths        : int [1:71] 11 11 1 2 1 2 1 2 1 3 ...
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  ..@ assayData        :<environment: 0x104fe810> 
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
