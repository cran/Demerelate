% File src/library/base/man/weir.Rd
     \name{weir}
     \alias{weir}
     \title{Calculation of \eqn{F_{is}} based on Weir and Cockerham 1984}
     \description{
       Function used internally by \code{Demerelate} and \code{F.stat} to calculate \eqn{F_{is}} values based on Weir and Cockerham (1984) statistics.
     }
     \usage{
     weir(tab.pop, tab.freq, popsize)
     }
     \arguments{
       \item{tab.pop}{Dataframe following format of \link{inputformat} only column three and four are used for calculations.}
       \item{tab.freq}{Vector of allele frequencies with name of each allele in column name.}
       \item{popsize}{Population size.}
		}

      \details{The function provides the alternative methods to Nei's \eqn{F_{is}} calculation. \eqn{F_{is}} values are calculated according to Weir and Cockerham 1984. Consult references for details}

      \value{
	\item{fis.weir}{\eqn{F_{is}} value calculated by method of Weir and Cockerham 1984}
}
	

\references{
    Weir, B.S. and Cockerham, C.C. (1984) Estimating F-Statistics for the analysis of population structure. Evolution, 38, 1358-1370.
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{Fis.calc}
	     \link{F.stat}
	     \link{Fis}
	     }
     }
     \examples{
     ## internal function of F.stat not intended for direct usage
     
     }
