% File src/library/base/man/allele.sharing.Rd
     \name{allele.sharing}
     \alias{allele.sharing}
     \title{Calculates allele sharing rates or similarity estimators for two populations}
     \description{
       Internal function of Demerelate to calculate different allele sharing indices or estimators for any pair of population. 
     }
     \usage{
     allele.sharing(pop1, pop2, allele.column, data,
                    value = NA, ref.pop)
     }
     \arguments{
       \item{pop1}{Specific dataframe of type \link{inputformat}. Population one used for calculations. Inputformat needs to be standard three digits mode for Demerelate.}
       \item{pop2}{Specific dataframe of type \link{inputformat}. Population two used for calculations. Inputformat needs to be standard three digits mode for Demerelate.}
       \item{allele.column}{Numeric value - It equals the number of the first column in the dataframe containing allele information. Order of loci in both populations needs to be exactly equal.}
       \item{data}{If set as TRUE allele sharing is only calculated for each pairwise comparison of individuals in the same position in pop1 with this position in pop2 - i.e. the diagonal of a distance matrix. If set as data.frame each pairwise comparison of column 1 with column 2 is calculated - i.e. the lower triangle of a distance matrix.}
       \item{value}{Character defining method to calculate allele sharing or similarity estimates. Can be set as "Bxy", "Sxy", "Mxy", "Li", "rxy", "lxy", "loiselle", "wang.fin", "wang", "ritland", "morans.fin" or "morans".}
       \item{ref.pop}{Reference population used for relatedness calculations.}
     }

      \details{The function uses each population to build a dataframe with row positions for pairwise comparisons according to their position in \code{pop1} and \code{pop2}. Similarity indices are calculated for each pair in one locus defined by the number indicated in \code{allele.column}. If \code{allele.column} is set as \code{1} allele information will be taken from the two columns of the first locus.}

      \value{Function returns the object simil.vector as vector of similarities according to \code{value}.}

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{inputformat}
	     \link{Demerelate}
	     \link{Emp.calc}
	     \link{Bxy}
	     \link{Sxy}
	     \link{Mxy}
	     \link{Li}
	     \link{rxy}
	     \link{lxy}
	     \link{loiselle}
	     \link{wang.fin}
	     \link{wang}
	     \link{ritland}
	     \link{morans.fin}
	      }
     }
     \examples{
     ## internal function not intended for direct use
     
     }
