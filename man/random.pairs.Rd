% File src/library/base/man/random.pairs.Rd
     \name{random.pairs}
     \alias{random.pairs}
     \title{Random pairs generator}
     \description{
       Generating of random pairs from a single table with format \link{inputformat}.
     }
     \usage{
     random.pairs(tab.all, allele.column, pairs)
     }
     \arguments{
		  \item{tab.all}{Dataframe containing genetic information with format \link{inputformat}.} 
		  \item{allele.column}{Loci information in column \code{allele.column} and \code{allele.column}+1.} 
		  \item{pairs}{Number of pairs to be generated.} 
     }

      \value{
		  \item{rand.pairs}{List of random pairs.}
		  
	      }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{offspring}
	      }
     }
     \examples{
     ## internal function not intended for direct usage 
    
     
     }
