% File src/library/base/man/demerelpop.Rd
     \name{demerelpop}
     \docType{data}
     \alias{demerelpop}
     \alias{inputformat}
     \title{Example dataset for package Demerelate to calculate inter individual pairwise genetic relatedness.}
     \description{
       The dataset gives randomized loci information of populations with differently related individuals. \emph{Pop-FS-HS} consists of 10
	full siblings (45 pairs) combined with 10 half siblings (45 pairs). \emph{Pop-FS-Non} consists of 10 full sibilings (45 pairs) and 10 random individuals (45 pairs). \emph{Pop-Non} is a population of 20 random drawn individuals (190 pairs). Datasets are based on information of 8 diploid loci with total number in alleles indicated by column names, i.e. number of alleles = 5, 10, 15, 20, 25, 30, 35, 40.
     }

     \usage{
     data(demerelpop)
     }

\details{ Standard inputformat is given for all applications in a similar way. Table is formatted as dataframe with headers set as TRUE. Headernames are not necessary but recommended.\cr 
First column contains sample IDs (mode=character), which should be unique. However, no proofreading is implemented yet, but double named individuals may lead to errors or strange results. Column two contains population information (mode=factor). 
Pairwise relatedness is only calculated within populations. If you want to compare pairwise individual comparisons of several populations you need to build a new fictive population. You should do that with caution, since hidden population structures may produce strange results. Column three and four containing two alleles of a diploid marker for the sample defined in column one. Each marker should be appended in pairwise columns from then on.\cr Allele size should be sorted in columns with the small allele in the odd numbered column and the bigger in even numbered column for convenience. However, Demerelate will take care of not sorted allele sizes and handle it correctly. The same inputformat is used for reference populations.\cr
All populations of inputdata are used as reference if no reference is defined by \code{reference.pop}.
Additionally, the distance data are given (\code{tab.dist}) by the same inputformat. In column three and four you will need to define the x and y coordinate of each individual either by relative x-y coordinates or by decimal degrees of geographic coordinates.

Example input. Every dataframe in Demerelate should be organized in this way.
\tabular{lllllll}{
Sample-ID\tab Population \tab locus.1.a \tab locus.1.b \tab locus.2.a \tab locus.2.b \tab ...\cr
Ind.Norway.01\tab Norway \tab 001 \tab 002 \tab 001 \tab 002 \tab ...\cr
Ind.Norway.02\tab Norway \tab 001 \tab 003 \tab 002 \tab 005 \tab ...\cr
Ind.Norway.03\tab Norway \tab 001 \tab 004 \tab 003 \tab 004 \tab ...\cr
Ind.Norway.04\tab Norway \tab 002 \tab 005 \tab 006 \tab 008 \tab ...\cr
Ind.Germany.01\tab Germany \tab 001 \tab 001 \tab 001 \tab 006 \tab ...\cr
Ind.Germany.02\tab Germany \tab 002 \tab 002 \tab 001 \tab 007 \tab ...\cr
Ind.Germany.03\tab Germany \tab 001 \tab 006 \tab 001 \tab 004 \tab ...\cr
Ind.Germany.04\tab Germany \tab 003 \tab 004 \tab 001 \tab 002 \tab ...\cr
}
}
     

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{Demerelate}
	     \link{demereldist}
	     \link{demerelref}
	      }
     }
     \examples{

# Please consult examples from Demerelate   
     
     }
