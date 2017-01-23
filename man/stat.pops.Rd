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
               pm, genotype.ref)
     }
     \arguments{
       \item{Thresholds}{Thresholds of relatedness to be used for statistics on empirical data.}
       \item{tab.pop.pop}{Dataframe following format of \link{inputformat} only column three and four are used for calculations.}
       \item{pairs}{Number of randomized pairings used for reference populations.}
       \item{p.correct}{\emph{logical} - should Yates correcture be used for \eqn{\chi^2} statistics.}
\item{directory.name}{Name of the directory results send to.} 
\item{out.name}{Filename of the output.}
\item{file.output}{\emph{logical}. Should a cluster dendogram, histograms and .txt files be sent as standard output in your working directory. In some cases (inflating NA values) it may be necessary that this value has to be set as FALSE due to problems in calculating clusters on pairwise NA values.} 
\item{inputdata}{R object or external file to be read internally with standard Demerelate \link{inputformat}. Dataframe will be split by population information and calculations will run separately. If no reference population information is specified (\code{ref.pop} = "NA") all information on loci are used as reference by omitting population information.}
\item{object}{Information whether inputdata are an object or should be read in as file.}
\item{value}{String defining method to calculate allele sharing or similarity estimates.}
\item{iteration}{Number of bootstrap iterations in \eqn{F_{is}} calculations.} 
\item{pm}{R object or external file to be read internally. Custom reference populations will be loaded for the analysis. The object should be a list of allele frequencies. For each locus a vector with allele frequencies p and allele names as vector names needs to be combined to a list. The last list object is a vector of sample sizes for each locus.}
\item{genotype.ref}{logical. If set as TRUE random non related populations are generated from genotypes of the reference population. If set as false allele frequencies are used for reference population generation.}

		}

      \details{Values for logical operators and data are given by \link{Demerelate}.}

      \value{
	\item{ClusterPopulationNameOutName.pdf}{Combined information on dataset by different types of plots.}
	\item{Relate.meanPopulationNameOutName.txt}{Combined information on dataset regarding calculated thresholds and number of different types of relatives in population. Summary of \eqn{\chi^2} statistics.}
	\item{out.stat[[1]]}{Empirical relatedness in populations.}
	\item{out.stat[[2]]}{Summary statistics for \eqn{\chi^2} statistics.}
  \item{error.individuals.txt}{If too many NAs are found in the data the program may stop here and tells you what individuals caused the troubles in calculations.}
}

	

\references{
     Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
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
