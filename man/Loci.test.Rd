% File src/library/base/man/Loci.test.Rd
     \name{Loci.test}
     \alias{Loci.test}
     \title{Analysis on differences in mean relatedness based on number of loci used for calculations}
     \description{
       Random samples are drawn from populations specified in \code{tab.pop}. Pairwise relatedness is calculated by either \emph{rxy}, \emph{Mxy} or \emph{Bxy}. Pairwise relatedness is calculated for different number of loci beginning with only one up to the maximum number of loci in the dataset. Reference populations can be set for the calculations if \code{ref.pop="NA"}. \code{tab.pop} is used as an estimate of overall allele frequencies. Distance between pairwise relatedness estimates are calculated using \code{dist()} with euclidean calculation.\cr Graphically, a .pdf file is exported and the object 
\code{Random.loci.differences} contains each pairwise difference for usage in further statistics. 
     
       
     }
     \usage{
     Loci.test(tab.pop, ref.pop = "NA", object = FALSE,
               value = "rxy", bt = 1000, 
               file.output = FALSE)
     }
     \arguments{
       \item{tab.pop}{Specific dataframe of type \link{inputformat}. All populations which should be analyzed.}
       \item{ref.pop}{Specific dataframe of type \link{inputformat}. Population information which should be used as a reference for either threshold or allele frequency calculation.}
       \item{object}{\emph{logical} - is tab.pop object.}
       \item{value}{String defining method to calculate allele sharing or similarity estimates. Can be set as "rxy", "Bxy" or "Mxy".}
       \item{bt}{Number of bootstrap replicates.}
\item{file.output}{\emph{logical} - should a .pdf file be sent as standard output in your working directory.}
            }

      \value{		
		  \item{Loci.test.Sys.Date().pdf}{Plot of linear regression of used loci and calculated mean estimate of relatedness from each bootstrapped replicate defined by \code{bt}.}
		  \item{Random.loci.differences}{List of differences in relatedness derived from bootstrap replicates. The x'th list element contains a matrix of means from replicates each over x random loci.}

	      }

\references{
    Blouin, M.S. et al. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.\cr
    Li C.C. and Horvitz D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics 5, 107-17.\cr
    Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
    Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.\cr
    Maechler, M. and many others. (2012) Utilities from Seminar fuer Statistik ETH Zurich. R package version 1.0-20.
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{inputformat}
	     \link{allele.sharing}
	     }
     }
     \examples{
     ## Loci.test to analyse Loci distribution in example data
     ## Bxy, rxy and Mxy are possible as estimators.

     data(demrelpop)
     demrelpop.sp <- split(demrelpop,demrelpop[,2])
     
     Loci.results <- Loci.test(demrelpop.sp[[1]][,1:8], 
                     ref.pop = "NA", object = TRUE,
                     value = "rxy", bt = 10)  
     
     }
