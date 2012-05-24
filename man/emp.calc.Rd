% File src/library/base/man/emp.calc.Rd
     \name{emp.calc}
     \alias{emp.calc}
     \title{Function to calculate pairwise relatedness within populations based on allele sharing.}
     \description{
       Allele sharing and relatedness estimators can be calculated for one population as a mean over several loci.
     }
     \usage{
     emp.calc(tab.pop.pop, value="rxy", ref.pop="NA")
     }
     \arguments{
       \item{tab.pop.pop}{Object with information for one population and loci. Object needs to be formatted according to \link{inputformat}}
       \item{value}{Either "rxy","Bxy","Mxy" are available to calculate mean relatedness among loci using \link{allele.sharing}}
\item{ref.pop}{Reference population needs to be specified for rxy calculations}
     }

      \details{The function calculates pairwise relatedness for all individuals in the dataframe. The output is a matrix of similarities by relatedness values in the population.}

      \value{
		  \item{empirical.list}{Object containing matrix of mean individual pairwise relatedness.}   

	      }

\references{
    Blouin, M.S. et al. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.\cr
    Li C.C. and Horvitz D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics 5, 107-17.\cr
    Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
    Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{allele.sharing}
	     \link{queller}
	      }
     }
     \examples{
     ## demrelpop needs to be an object of class data.frame formatted
     ## according to inputformat. For each pairwise comparison in
     ## in the data.frame empirical estimates of relatedness are
     ## calculated omitting statistics.
     ## Bxy, rxy and Mxy are possible to chose as estimator.

      data(demrelpop)
      
      # As example Mxy is calculated for one population of demrelpop
      
      demrelpop.sp <- split(demrelpop,demrelpop[,2])

      empirical.result <- emp.calc(demrelpop.sp[[1]], value="Mxy", 
                                   ref.pop="NA")
    
     
     }
