% File src/library/base/man/Demerelate.Rd
     \name{Demerelate}
     \alias{Demerelate}
     \title{Demerelate --- Algorithms to estimate pairwise relatedness within populations based on allele sharing}
     \description{
       Head function of Demerelate. This function should be called if any estimation of relatedness is intended. Additionally, some F-statistics can be calculated. Default parameters are set for convenient usage. Only an input dataframe containing allelic information is necessary. Geographical distances, reference populations or alterations on statistics can be set by adapting parameters.}
       
       
     \usage{
     Demerelate(inputdata, tab.dist = "NA", ref.pop = "NA", 
                object = FALSE, value = "Mxy", Fis = FALSE,
                file.output = FALSE, p.correct = FALSE,
                iteration = 1000, pairs = 1000, 
                dis.data = "relative", NA.rm = TRUE,
                genotype.ref = TRUE)
     }
     
     \arguments{
       \item{inputdata}{R object or external file to be read internally with standard Demerelate \link{inputformat}. Dataframe will be split by population information and calculations will run separately. If no reference population information is specified (\code{ref.pop = "NA"}) all information on loci are used as reference by omitting population information.}
       \item{tab.dist}{R object or external file to be read internally with standard Demerelate \link{inputformat}. Geographic distances can be defined and will be analysed combined with genetic data. Column three and four of standard inputformat are used for x and y coordinates.}
       \item{ref.pop}{R object or external file to be read internally with standard Demerelate \link{inputformat}. Custom reference populations will be loaded for the analysis. Population information of reference file will be omitted so that allele frequencies are calculated from the whole dataset. Optionally allele frequencies can be loaded as reference: The object should be then a list of allele frequencies. For each locus a vector with allele frequencies p and allele names as vector names needs to be combined to a list. The last list object is a vector of sample sizes for each locus.}
       \item{object}{Information whether inputdata are objects or should be read in as files.}
       \item{value}{String defining method to calculate allele sharing or similarity estimates. Can be
 set as "Bxy", "Sxy", "Mxy", "Li", "lxy", "rxy", "loiselle", "wang.fin", "wang",
  "ritland", "morans.fin" or "morans" \link{allele.sharing}.}
       \item{Fis}{\emph{logical}. Should \eqn{F_{is}} values be calculated for each population?}
       \item{iteration}{Number of bootstrap iterations in \eqn{F_{is}} calculations.}
       \item{pairs}{Number of pairs calculated from reference populations for randomized full siblings, half siblings and non related individuals.}
       \item{file.output}{\emph{logical}. Should a cluster dendogram, histograms and .txt files be sent as standard output in your working directory. In some cases (inflating NA values) it may be necessary that this value has to be set as FALSE due to problems in calculating clusters on pairwise NA values.}
       \item{p.correct}{\emph{logical}. Should Yates correction from \code{prop.test(...)} be used in \eqn{\chi^2} statistics when calculating p-values on differences between empirical and randomized relatedness in populations.}
       \item{dis.data}{The kind of data to be used as distance measure. Can be "relative" - relative x and y coordinates should be given in \code{tab.dist} or "decimal" for geographic decimal degrees.}
       \item{NA.rm}{logical. If set as TRUE samples with NA in any position are removed from the calculation. If set as FALSE you may get an error message telling you to remove some individuals to run through the procedure. Always be aware that if your calculations are successful although you have NA values in your populations your may be biased by missing data.}
       \item{genotype.ref}{logical. If set as TRUE random non related populations are generated from genotypes of the reference population. If set as false allele frequencies are used for reference population generation. If ref.pop is given as list of allele frequencies genotype.ref = FALSE is forced.}
     }

      \details{
      Pairwise relatedness is calculated from inputdata. Be sure to fit exactly the \link{inputformat}. Missing values are omitted when flagged as \code{NA}. If no additional reference populations are defined, inputdata omitting population information are used to calculate references. If no good reference populations are available you need to take care of bias in calculations. In any case you should consult for example Oliehoek et al. 2006 to get an idea of bias in relatedness calculations. \cr
      Geographic distances between individual pairs are calculated when \code{tab.dist =} ... . Distances calculated from x-y coordinates by simple Pythagorean mathematics can be applied to any metrical positions in sampling. Geographic coordinates from e.g. GPS need to be transformed to decimal GPS coordinates. Be sure to have positions for each individual or remove missing values from inputdata.\cr
      Each calculation will have its unique bar-code and is named with the date and population name. Calculations are performed for each population in the inputdata.
	}
      \value{
      
		\bold{Function returns files in a folder named with a bar-code and date of analysis as follows if file.output is set as TRUE:}
		  \item{Empirical.relatedness.Population.txt}{Matrix of relatedness values for each population.}
		  \item{Geographic.distances.Population.txt}{Matrix of geographic distances for each population.}
		  \item{Relate.mean.Populationout.name.txt}{Depends on selected estimators and mode of analysis. Either a summary of correlation of relatedness with geographic distance for each population or a summary of tests for relatedness within populations compared to reference populations is written to the file.}
		  \item{Random.Halfsib.distances.overall.txt}{Matrix of relatedness values calculated from randomized reference population for half siblings.}
		  \item{Random.NonRelated.distances.overall.txt}{Matrix of relatedness values calculated from randomized reference population for non related individuals.}
		  \item{Random.Fullsib.distances.overall.txt}{Matrix of relatedness values calculated from randomized reference population for full siblings.}
		  \item{Cluster.Populationout.name.pdf}{Containing an UPGMA cluster dendogram of relatedness values and a histogram of relatedness values per locus and for loci overall.}
		  \item{Total-Regression.Population.pdf}{Containing regression plot and linear fit for geographic distance and genetic relatedness.}
		  \item{Summary.Populationout.name.txt}{Summary of analysis of F statistics and allele/genotype frequencies.}
                  

      
		\bold{Function returns via \code{return} following objects as one list:}
                  
      \item{dem.results[[1]]}{Settings of the calculation are passed to this list object.}
		  \item{dem.results[[2]]}{Mean relatedness for empirical population over all loci.}
		  \item{dem.results[[3]]}{Summarized relatedness statistics with thresholds and randomized populations from the dataset.}
		  \item{dem.results[[4]]}{Statistical analysis of the number of siblings found for each population.}
 		  \item{dem.results[[5]]}{Thresholds for relatedness if "Bxy" or "Mxy" are selected as estimators}
		  \item{dem.results[[6]]}{\eqn{F_{is}} values and statistics for each population if \code{Fis==TRUE}}
		  \item{dem.results[[7]]}{Summary of linear regression of distance data are provided.}
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
1203-1215.\cr  }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{\link{inputformat}
	     \link{Emp.calc}
	     \link{stat.pops}
	     \link{F.stat}
		}
	     }

     \examples{   

     
     ## Data set is used to calculate Blouins allele sharing index on  
     ## population data. Pairs are set to 10 for convenience.
     ## For statistical reason for your final results you may want to 
     ## use more pairs to model relatedness (1000 pairs are recommended).

     data(demerelpop)
     
     dem.results <- Demerelate(demerelpop[,1:6], value="Mxy", 
                    file.output=FALSE, object=TRUE, pairs=10)


     ## Demerelate can be executed with several different values 
     ## should consult the references to decided which estimator may 
     ## be useful in your case. 
     ## Be careful some estimators may be biased in situations with
     ## no reference populations or violatin of Hardy-Weinberg
     ## Equilibrium.
    
     
     }
