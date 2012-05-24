% File src/library/base/man/Fis.calc.Rd
     \name{Fis.calc}
     \alias{Fis.calc}
     \title{Calculation of \eqn{Fis} empirical and bootstrapped values}
     \description{
       Internal function of Demerelate and F.stat to calculate \eqn{F_{is}} values from empirical datasets. Performing randomization statistics and preparing output.
     }
     \usage{
     Fis.calc(tab.pop, iteration, number.loci, object, 
              directory.name, out.name)
     }
     \arguments{
       \item{tab.pop}{Data.frame following format of \link{inputformat} only column three and four are used for calculations}
       \item{iteration}{Number of replicates for bootstrap statistics}
       \item{number.loci}{Number of loci in tab.pop}
       \item{object}{Whether tab.pop is an object or file}
       \item{directory.name}{Name of the directory results send to}
       \item{out.name}{Filename of the output}
       

		}

      \value{
	Function returns a list containing the following information:
       
	\item{output.fis[[1]]}{Empirical \eqn{F_{is}} values according to Nei 1977.}
	\item{output.fis[[2]]}{Empirical \eqn{F_{is}} values according to Weir and Cockerham 1984.}
	\item{output.fis[[3]]}{p values for significance for \eqn{F_{is}} values according to Nei 1977.}
        \item{output.fis[[4]]}{p values for significance for \eqn{F_{is}} values according to Weir and Cockerham 1984.}

Additionally, a file is generated containing all these information if file.output is set as TRUE or alternatively when directory.name and out.name are given.
 
	\item{SummaryPopulationout.name.txt}{Combined output with different \eqn{F_{is}} metrices and allele\\genotype frequencies} 
		}
	

\references{
    Weir, B.S. and Cockerham, C.C. (1984) Estimating F-Statistics for the analysis of population structure. Evolution, 38, 1358-1370.\cr
    Nei, M. (1977) F-statistics and analysis of gene diversity in subdivided populations. Annals of Human Genetics, 41, 225-233.
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{weir}
	     \link{F.stat}
	     \link{Fis}
	     }
     }
     \examples{
     ## internal function of F.stat not intended for direct usage
     
     }
