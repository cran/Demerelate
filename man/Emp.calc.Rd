% File src/library/base/man/Emp.calc.Rd
     \name{Emp.calc}
     \alias{Emp.calc}
     \title{Function to calculate pairwise relatedness within populations based on allele sharing.}
     \description{
       Allele sharing and relatedness estimators can be calculated for one population as a mean over several loci.
     }
     \usage{
     Emp.calc(tab.pop.pop, value="NA", ref.pop="NA")
     }
     \arguments{
       \item{tab.pop.pop}{Object with information for one population and loci. Object needs to be formatted according to \link{inputformat}.}
       \item{value}{String defining method to calculate allele sharing or similarity estimates. Can be
 set as "Bxy", "Sxy", "Mxy", "Li", "lxy", "rxy", "loiselle", "wang.fin", "wang",
  "ritland", "morans.fin" or "morans" \link{allele.sharing}.}
\item{ref.pop}{Reference population needs to be specified for several calculations.
}
     }

      \details{The function calculates pairwise relatedness for all individuals in the dataframe. The output is a matrix of similarities by relatedness values in the population.}

      \value{
		  \item{empirical.list}{Object containing a dataframe of mean individual pairwise relatedness. If value="NA" all available estimators are calculated in one call of function \code{Emp.calc}.}   

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
       \code{
	     \link{allele.sharing}
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
	     \link{morans}
	      }
     }
     \examples{
     ## demerelpop needs to be an object of class data.frame formatted
     ## according to inputformat. For each pairwise comparison in
     ## in the data.frame empirical estimates of relatedness are
     ## calculated omitting statistics.

      data(demerelpop)
      
      # As example Mxy is calculated for one population of demerelpop
      
      demerelpop.sp <- split(demerelpop,demerelpop[,2])

      empirical.result <- Emp.calc(demerelpop.sp[[1]], value="Mxy", 
                                   ref.pop="NA")
    
     
     }
