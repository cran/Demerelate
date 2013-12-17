% File src/library/base/man/allele.sharing.Rd
     \name{allele.sharing}
     \alias{allele.sharing}
     \title{Calculates allele sharing rates or similarity estimators for two populations}
     \description{
       Internal function of Demerelate to calculate different allele sharing indices or estimators for any pair of population. 
       Following indices can be calculated: 
	    \eqn{B_{xy}} (number of alleles shared) as described in Li and Horvitz 1953. 
	    \eqn{M_{xy}} (genotype sharing) as described in Blouin et al. 1996. 
	    The estimator \eqn{r_{xy}} based on Queller and Goodnight 1989 adapted to pairwise comparisons as described in Oliehoek et al. 2006.
     }
     \usage{
     allele.sharing(pop1, pop2, allele.column, onlypairs = FALSE,
                    value = NA, ref.pop)
     }
     \arguments{
       \item{pop1}{Specific dataframe of type \link{inputformat}. Population one used for calculations. Individuals passed to rows of resulting matrix. Inputformat needs to be standard three digits mode for Demerelate.}
       \item{pop2}{Specific dataframe of type \link{inputformat}. Population two used for calculations. Individuals passed to columns of resulting matrix. Inputformat needs to be standard three digits mode for Demerelate.}
       \item{allele.column}{Numeric value which needs to be an odd number. It equals the number of the first column in the dataframe containing allele information. Order of loci in both populations needs to be exactly equal.}
       \item{onlypairs}{If set as TRUE allele sharing is only calculated for diagonal comparisons in the resulting pairwise relatedness matrix. If set as FALSE lower triangle is calculated from the pairwise relatedness matrix.}
       \item{value}{Character defining method to calculate allele sharing or similarity estimates. Can be set as "rxy", "Bxy" or "Mxy".}
       \item{ref.pop}{Reference population used for relatedness calculations.}
     }

      \details{The function uses column one of each population dataframe to build a matrix with column names from samples of \code{pop2} and row names from samples of \code{pop1}. Similarity indices can be calculated for each pair. Indices are calculated using allelic information of one locus which first column must be indicated by \code{allele.column}. \code{allele.column}+1 will be used automatically as second allele information.}

      \value{Function returns the object tab.all as matrix of rows passed from population one and columns from population two.}

\references{
    Blouin, M.S. et al. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.\cr
    Li C.C. and Horvitz D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics 5, 107-17.\cr
    Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
    Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{inputformat}
	     \link{Demerelate}
	     \link{emp.calc}
	     \link{queller}
	      }
     }
     \examples{
     ## internal function not intended for direct use
     
     }
