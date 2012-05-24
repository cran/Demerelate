% File src/library/base/man/offspring.Rd
     \name{offspring}
     \alias{offspring}
     \title{Mendelian offpsring generator}
     \description{
       Random generation of offspring from two parent individuals for one locus.
     }
     \usage{
     offspring(parent1, parent2, allele.column, pairs)
     }
     \arguments{
		  \item{parent1}{data.frame of parent one.} 
		  \item{parent2}{data.frame of parent two.}
		  \item{allele.column}{Loci information in column \code{allele.column} and \code{allele.column}+1.} 
		  \item{pairs}{Number of offspring generated from parent1 and parent2.} 
     }

      \details{Mendelian generator of random offspring.}

      \value{
		  \item{offspr}{Dataframe of offpsring generated with length pairs.}
		  
	      }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{random.pairs}
	      }
     }
     \examples{
     ## internal function not intended for direct usage
    
     
     }
