% File src/library/base/man/Loci.test.Rd
     \name{Loci.test}
     \alias{Loci.test}
     \title{Analysis on differences in mean relatedness based on number of loci used for calculations}
     \description{
       Random samples are drawn from populations specified in \code{tab.pop}. Pairwise relatedness is calculated by either "Bxy", "Sxy", "Mxy", "Li", "lxy", "rxy", "loiselle", "wang.fin", "wang",
  "ritland", "morans.fin" or "morans" \link{allele.sharing}. Pairwise relatedness is calculated for different number of loci beginning with only one up to the maximum number of loci in the dataset. Reference populations can be set for the calculations if \code{ref.pop="NA"}. \code{tab.pop} is used as an estimate of overall allele frequencies. Distance between pairwise relatedness estimates are calculated using \code{dist()} with euclidean calculation.\cr Graphically, a .pdf file is exported and the object \code{Random.loci.differences} contains each pairwise difference for usage in further statistics. 
     
       
     }
     \usage{
     Loci.test(tab.pop, ref.pop = "NA", object = FALSE,
               value = NA, bt = 1000, 
               file.output = FALSE)
     }
     \arguments{
       \item{tab.pop}{Specific dataframe of type \link{inputformat}. All populations which should be analyzed.}
       \item{ref.pop}{Specific dataframe of type \link{inputformat}. Population information which should be used as a reference for either threshold or allele frequency calculation.}
       \item{object}{\emph{logical} - is tab.pop object.}
       \item{value}{String defining method to calculate allele sharing or similarity estimates. Can be set as "Bxy", "Sxy", "Mxy", "Li", "lxy", "rxy", "loiselle", "wang.fin", "wang", "ritland", "morans.fin" or "morans".
}
       \item{bt}{Number of bootstrap replicates.}
\item{file.output}{\emph{logical} - should a .pdf file be sent as standard output in your working directory.}
            }

      \value{		
		  \item{Loci.test.Sys.Date().pdf}{Plot of linear regression of used loci and calculated mean estimate of relatedness from each bootstrapped replicate defined by \code{bt}.}
		  \item{Random.loci.differences}{List of differences in relatedness derived from bootstrap replicates. The x'th list element contains a matrix of means from replicates each over x random loci.}

	      }

\references{
        Armstrong, W. (2012) fts: R interface to tslib (a time series library in c++). by R package version
0.7.7.\cr
Blouin, M., Parsons, M., Lacaille, V. and Lotz, S. (1996) Use of microsatellite loci to classify indi-
viduals by relatedness. Molecular Ecology, 5, 393-401.\cr
Hardy, O.J. and Vekemans, X. (1999) Isolation by distance in a contiuous population: reconciliation
between spatial autocorrelation analysis and population genetics models. Heredity, 83, 145-154.\cr
Li, C.C., Weeks, D.E. and Chakravarti, A. (1993) Similarity of DNA fingerprints due to chance and
relatedness. Human Heredity, 43, 45-52.\cr
Li, C.C. and Horvitz, D.G. (1953) Some methods of estimating the inbreeding coefficient. Ameri-
can Journal of Human Genetics, 5, 107-17.\cr
Loiselle, B.A., Sork, V.L., Nason, J. and Graham, C. (1995) Spatial genetic structure of a tropical
understory shrub, Psychotria officinalis (Rubiaceae). American Journal of Botany, 82, 1420-1425.\cr
Lynch, M. (1988) Estimation of relatedness by DNA fingerprinting. Molecular Biology and Evolu-
tion, 5(5), 584-599.\cr
Lynch, M. and Ritland, K. (1999) Estimation of pairwise relatedness with molecular markers. Ge-
netics, 152, 1753-1766.\cr
Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with
a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution,
43, 258-275.\cr
Ritland, K. (1999) Estimators for pairwise relatedness and individual inbreeding coefficients. Ge-
netics Research, 67, 175-185.\cr
Wang, J. (2002) An estimator for pairwise relatedness using molecular markers. Genetics, 160,
1203-1215.\cr 
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{inputformat}
	     \link{allele.sharing}
	     }
     }
     \examples{
     ## Loci.test to analyse Loci distribution in example data

     data(demerelpop)
     demerelpop.sp <- split(demerelpop,demerelpop[,2])
     
     Loci.results <- Loci.test(demerelpop.sp[[1]][,1:8], 
                     ref.pop = "NA", object = TRUE,
                     value = "rxy", bt = 10)  
     
     }
