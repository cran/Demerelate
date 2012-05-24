% File src/library/base/man/Lin.reg.distance.Rd
     \name{Lin.reg.distance}
     \alias{Lin.reg.distance}
     \title{Linear regression of empirical genetic relatedness with geographic distance.}
     \description{
       Function caluclates CI intervals and exports plots and matrices for linear regression of empirical genetic relatedness with geographic distance.
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

      \details{Regression model is calculated by lm(...) for datasets. Be aware that your data need to fit assumptions for linear regression statistics. You should consult more statistics to check for e.g. normality or heteroskedastics in your data. Linear regression may be in most cases useful if calculating \eqn{r_{xy}} as estimator. The \eqn{B_{xy}} and \eqn{M_{xy}} may need considerable sample sizes to fit general assumptions for tests.}

      \value{
		  \item{Total-Regression.Population.pdf}{Containing regression plot and linear fit for geographic distance and genetic relatedness.}
		  \item{Relate.mean.Populationout.name.txt}{A summary of correlation of relatedness with geographic distance.}
		  
		  Summary statistics from linear regression of geogrphical distance and mean relatedness over loci are exported as list object.
		  
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
