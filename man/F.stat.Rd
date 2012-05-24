% File src/library/base/man/F.stat.Rd
     \name{F.stat}
     \alias{F.stat}
     \title{Head function to calculate F statistics.}
     \description{Calculation of \eqn{F_{is}} values for a single population and several loci. Calculations are based on statistics from Weir and Cockerham 1984 and Nei 1977.}

     \usage{
     F.stat(tab.pop, object = TRUE, iteration = 1000, 
            directory.name = "NA", out.name = "NA")
     }

     \arguments{
		  \item{tab.pop}{A file with a dataframe containing genetic information with format \link{inputformat}.}
		  \item{object}{logical - whether data are object or file.}  
		  \item{iteration}{Number of bootstrap replicates} 
		  \item{directory.name}{Directory name where files should be exported. Directory must be present when running the function. If set as "NA" results are only printed on screen.} 
		  \item{out.name}{Name of file for exporting summary of statistics.} 
     }

      \details{This function is normally executed by Demerelate when \eqn{F_{is}} set as TRUE. However, \eqn{F_{is}} values can be calculated on single populations without the need of doing a whole new relatedness analysis. This may be useful if \eqn{F_{is}} values were not calculated in the first place. Be sure that the path already is present prior to analysis. However, if set as NA results are printed on screen or exported to an object. Inputdata will be splitted by population names and calculations will be made for each single population. Note that mean \eqn{F_{is}} values calculated by Weir and Cockerham's method are either given as arithmetic mean of single loci \eqn{F_{is}} and as weighted average over loci, which should be unbiased to a first aproximation according to Weir and Cockerham 1984. Since there are many ways to calculate mean \eqn{F_{is}} and all may serve in different situations as 'good' estimate you should spent some thoughts on this issue.}

      \value{
      
        Function returns a list for each population containing the following information:
       
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
       \code{
	     \link{Fis.calc}
	      }
     }

     \examples{

        ## Fis statistics are calculated on demrelpop. Weir and Cockerham
        ## and Nei estimates are provided and either written to an object
        ## or to a file.
        ## For end results iteration of 1000 is recommended.


	      data(demrelpop)
	      fstat.results <- F.stat(demrelpop, iteration = 10, 
                                directory.name = "NA", 
                                out.name = "NA")
        
     }
