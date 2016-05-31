% File src/library/base/man/relate.calc.Rd
     \name{relate.calc}
     \alias{relate.calc}
     \title{Coordinates internally reference and empirical datasets and statistics}
     \description{
       Internal function of Demerelate to combine different populations of randomized offspring and empirical populations to calculate thresholds and statistics for later use as reference for each empirical population.
     }
     \usage{
     relate.calc(tab.pop, pairs, file.output, value, 
                 directory.name, ref.pop)
     }
     \arguments{
       \item{tab.pop}{Dataframe following format of \link{inputformat} used as reference for randomizations.}
       \item{pairs}{Number of pairwise comparison for each randomization.}
\item{file.output}{\emph{logical}. Should a cluster dendogram, histograms and .txt files be sent as standard output in your working directory. In some cases (inflating NA values) it may be necessary that this value has to be set as FALSE due to problems in calculating clusters on pairwise NA values.} 
\item{value}{String defining method to calculate allele sharing or similarity estimates. Can be set as "Bxy", "Sxy", "Mxy", "Li", "lxy", "rxy", "loiselle", "wang.fin", "wang", "ritland", "morans.fin" or "morans".}
\item{directory.name}{Name of the directory results send to.}
\item{ref.pop}{R object or external file to be read internally with standard \code{Demerelate} \link{inputformat}. Custom reference populations will be loaded for the analysis. Population information of reference file will be omitted so that allele frequencies are calculated from the whole dataset.}
		}

      \details{The function internally calls any procedure of randomization and offspring generation. Finally, it coordinates the preparation of multinomial logistic regression for threshold calculation.}

      \value{
	\item{relate.return[[1]]}{\code{Object} - matrix of mean pairwise relatedness of all loci of each comparison in tab.pop.}
	\item{relate.return[[2]]}{\code{Object} - matrix of mean pairwise relatedness of each randomized comparison of full siblings.}
	\item{Random.Fullsib.relatedness.overall.txt}{\code{file} - matrix of mean pairwise relatedness of each randomized comparison of full siblings.}
	\item{relate.return[[3]]}{\code{Object} - matrix of mean pairwise relatedness of each randomized comparison of half siblings.}
	\item{Random.Halfsib.relatedness.overall.txt}{\code{file} - matrix of mean pairwise relatedness of each randomized comparison of half siblings.}
	\item{relate.return[[4]]}{\code{Object} - matrix of mean pairwise relatedness of each randomized comparison of non relatives.}
	\item{Random.NonRelated.relatedness.overall.txt}{\code{file} - matrix of mean pairwise relatedness of each randomized comparison of non relatives.}
	\item{relate.return[[5]]}{\code{Object} - calculated thresholds for relatedness from reference populations.}

}

\references{
    Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{Emp.calc}
	     \link{offspring}
	     \link{random.pairs}
	     \link{glm.prep}
	      }
     }
     \examples{
     ## internal function not intended for usage without Demerelate
     
     }
