% File src/library/base/man/stat.pops.Rd
     \name{stat.pops}
     \alias{stat.pops}
     \title{Calculation of Fis empirical and bootstrapped values}
     \description{
       Internal function of Demerelate to use prepared thresholds of relatedness to calculate state of relatedness for empirical data. Additionally, it combines information in different plots and prepares each single output.
     }
     \usage{
     stat.pops(Thresholds, tab.pop.pop, pairs, p.correct, 
               directory.name, out.name, file.output, 
               inputdata, object, value, iteration, 
               ref.pop)
     }
     \arguments{
       \item{Thresholds}{Thresholds of relatedness to be used for statistics on empirical data.}
       \item{tab.pop.pop}{Dataframe following format of \link{inputformat} only column three and four are used for calculations.}
       \item{pairs}{Number of randomized pairings used for reference populations.}
       \item{p.correct}{\emph{logical} - should Yates correcture be used for \eqn{\chi^2} statistics.}
\item{directory.name}{Name of the directory results send to.} 
\item{out.name}{Filename of the output.}
\item{file.output}{\emph{logical}. Should a cluster dendogram, histograms and .txt files be sent as standard output in your working directory. In some cases (inflating NA values) it may be necessary that this value has to be set as FALSE due to problems in calculating clusters on pairwise NA values.} 
\item{inputdata}{R object or external file to be read internally with standard Demerelate \link{inputformat}. Dataframe will be split by population information and calculations will run separately. If no reference population information is specified (reference.pop = "NA") all information on loci are used as reference by omitting population information.}
\item{object}{Information whether inputdata are an object or should be read in as file.}
\item{value}{String defining method to calculate allele sharing or similarity estimates. Can be set as "rxy", "Bxy" or "Mxy".}
\item{iteration}{Number of bootstrap iterations in \eqn{F_{is}} calculations.} 
\item{ref.pop}{R object or external file to be read internally with standard Demerelate \link{inputformat}. Custom reference populations will be loaded for the analysis. Population information of reference file will be omitted so that allele frequencies are calculated from the whole dataset.}

		}

      \details{Values for logical operators and data are given by \link{Demerelate}.}

      \value{
	\item{ClusterPopulationNameOutName.pdf}{Combined information on dataset by different types of plots.}
	\item{Relate.meanPopulationNameOutName.txt}{Combined information on dataset regarding calculated thresholds and number of different types of relatives in population. Summary of \eqn{\chi^2} statistics.}
	\item{out.stat[[1]]}{Empirical relatedness in populations.}
	\item{out.stat[[2]]}{Summary statistics for \eqn{\chi^2} statistics.}
}

	

\references{
     Blouin, M.S. et al. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.\cr
     Li C.C. and Horvitz D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics 5, 107-17.\cr
     Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
     Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.
	    }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{weir}
	     \link{F.stat}
	     }
     }
     \examples{
     ## internal function not intended for usage without Demerelate
     
     }
