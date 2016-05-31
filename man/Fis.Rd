% File src/library/base/man/Fis.Rd
     \name{Fis}
     \alias{Fis}
     \title{Calculates allele and genotype frequencies}
     \description{Internal function to prepare allele and genotype frequencies for F statistics. Rows with \code{NA} values are removed from the calculation.}

     \usage{Fis(tab.pop, allele.column)}

     \arguments{
       \item{tab.pop}{Object with information for one or several populations and loci. Format needs to be in accordance with \link{inputformat}.}
       \item{allele.column}{Loci information from column \code{allele.column} and \code{allele.column}+1 are used for the calculations.}
     }

      \details{Be sure to remove all missing data in your data file. Prior to calculation of \eqn{F_{is}} values. Every sample containing \code{NA} values in any column will be removed. It may be better to remove complete loci columns to get good results. In any case you should check this.}

      \value{
		  \item{tab.freq}{Table of allele frequencies}
		  \item{tab.freq.gen}{Table of genotype frequencies}
		  \item{Homo}{Sum of homozygotes}
		  \item{pop.size}{Population size}
		  \item{Fis}{\eqn{F_{is}} values calculated according to Nei 1983}
		  \item{Fis.weir}{\eqn{F_{is}} values calculated according to Weir and Cockerham 1984}
	      }

\references{
    Weir, B.S. and Cockerham, C.C. (1984) Estimating F-Statistics for the analysis of population structure. Evolution, 38, 1358-1370.\cr
    Nei, M. (1977) F-statistics and analysis of gene diversity in subdivided populations. Annals of Human Genetics, 41, 225-233.\cr
    Nei, M. and Chesser R.K. (1983) Estimation of fixation indices and gene diversities. Annals of Human Genetics, 47, 253-259.
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{weir}
	     \link{Fis.calc}
	     \link{F.stat}
	      }
     }
     \examples{
      ## internal function of F.stat not intended for direct use
     
     }
