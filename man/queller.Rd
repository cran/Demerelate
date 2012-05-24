% File src/library/base/man/queller.Rd
     \name{queller}
     \alias{queller}
     \title{Queller and Goodnight estimator of relatedness}
     \description{
       Calculation of symetric allele sharing using the RE-RAT method and allele frequencies for correction based on a reference population.
     }
     \usage{
     queller(pop1, pop2, allele.column, row.position,
             col.position, matrix.share, ref.pop)
     }
     \arguments{
		  \item{pop1}{Population data containing allele information. All individual names are found in row names of pairwise relatedness matrix.} 
		  \item{pop2}{Population data containing allele information. All individual names are found in column names of pairwise relatedness matrix.}
		  \item{allele.column}{Column number where loci information is taken from 
				       i.e. \code{allele.column} and \code{allele.column+1}.} 
		  \item{row.position}{Row position of the pairwise relatedness matrix which contains the name of one individual from population data \code{pop1}.} 
		  \item{col.position}{Column position of the pairwise relatedness matrix which contains the name of one individual from population data \code{pop2}.}
		  \item{matrix.share}{Pairwise relatedness matrix containing individual names of \code{pop1} and \code{pop2}.}
		  \item{ref.pop}{Reference population to calculate allele frequencies necessary for the calculation of the estimator.}
     }

      \details{Be aware that unlike in some other programs single individuals are not removed from the calculation to correct for any bias. If you want to use custom reference population to calculate overall allele frequencies, which may be corrected for single individuals or parts of the population, these population need to be constructed manually and indicated by ref.pop. Otherwise it may be feasible to use the whole population information to calculate allele frequency distributions.}

      \value{
		  \item{share}{value of \eqn{r_{xy}}. NA is returned if estimator could not be calculated.}
		  
	      }

\references{ 
    Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
    Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.}

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{random.pairs}
	      }
     }
     \examples{
     ## internal function not intended for direct usage 
    
     
     }
