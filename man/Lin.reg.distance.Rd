% File src/library/base/man/Lin.reg.distance.Rd
     \name{Lin.reg.distance}
     \alias{Lin.reg.distance}
     \title{Linear regression of empirical genetic relatedness with geographic distance.}
     \description{
       Function caluclates mantel statistics and exports plots and matrices for regression of empirical genetic relatedness with geographic distance.
     }
     \usage{
     Lin.reg.distance(dist.m, emp.dist, pairs, tab.pop.pop, 
                      offhalf.list, offfull.list, 
                      relate.off.non.Mxy.mean, file.output, 
                      directory.name, out.name, inputdata, 
                      object, value, iteration)
     }

     \arguments{
       \item{dist.m}{Object provided by \code{geo.dist}.}
       \item{emp.dist}{Empirical relatedness of population.}
       \item{pairs}{Number of bootstrap replicates to calculate confidence interval for linear regression.}
       \item{tab.pop.pop}{Population information to model linear fit.}
       \item{offhalf.list}{\code{Object} - matrix of mean pairwise relatedness of each randomized comparison of half siblings.}
       \item{offfull.list}{\code{Object} - matrix of mean pairwise relatedness of each randomized comparison of full siblings.}
       \item{relate.off.non.Mxy.mean}{\code{Object} - matrix of mean pairwise relatedness of each randomized comparison of non relatives..}
       \item{file.output}{\emph{logical} - Should a cluster dendogram, histograms and .txt files be sent as standard output in your working directory. In some cases (inflating NA values) it may be necessary that this value has to be set as FALSE due to problems in calculating clusters on pairwise NA values.}
       \item{directory.name}{Name of the directory results are send to.}
       \item{out.name}{Filename of the output.}
       \item{inputdata}{R object or external file to be read internally with standard Demerelate \link{inputformat}. Dataframe will be split by population information and calculations will run separately. If no reference population information is specified (reference.pop = "NA") all information on loci are used as reference by omitting population information.}
       \item{object}{Information whether inputdata are an object or should be read in as file.}
       \item{value}{String defining method to calculate allele sharing or similarity estimates. Can be set as "rxy", "Bxy" or "Mxy".}
       \item{iteration}{Number of bootstrap iterations in \eqn{F_{is}} calculations.}
     }

      \details{Mantel statistics is calculated by mantel(...) from vegan package for datasets. Similarities are transformed to distances by D=1-S for each estimator. For visualization a scatterplot is exported with relatedness thresholds calculated for each dataset to get an idea of the relatedness between genetic and geographic distance in addition to mantel statistics.}

      \value{
		  \item{Total-Regression.Population.pdf}{Containing regression plot and linear fit for geographic distance and genetic relatedness. Randomly generated populations for Full siblings, Halfsiblings and Non related individuals are used as optical reference in the plot. The mean relatedness of each population is plotted as line - blue dotted line for full siblings, red dotted line for half siblings and black dotted line for non related individuals.}
		  \item{Relate.mean.Populationout.name.txt}{A summary of correlation of relatedness with geographic distance.}
		  
		  Summary statistics from mantel correlation of geogrphical distance and mean relatedness over loci are exported as list object.
		  
	      }

\references{
      Oksanen, J. et al. (2013) vegan: Community Ecology Package. R package version 2.0-8
      }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{geo.dist}
	     \link{allele.sharing}
	     \link{Demerelate}
	      }
     }
     \examples{
    ## internal function not intended for direct use
    
     
     }
