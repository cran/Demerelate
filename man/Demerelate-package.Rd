\name{Demerelate-package}
\alias{Demerelate-package}
\docType{package}
\title{Demerelate --- Algorithms to estimate pairwise relatedness within populations based on allele sharing}

\description{
The package Demerelate provides several functions to calculate relatedness of individuals based on diploid genetic markers. Following indices can be calculated:

\enumerate{
\item \eqn{B_{xy}} (number of alleles shared) as described in Li and Horvitz 1953 (\link{Bxy}).
\item \eqn{S_{xy}} (number of alleles shared) as described in Lynch 1988 (\link{Sxy}).
\item  \eqn{M_{xy}} (genotype sharing) as described in Blouin et al. 1996 (\link{Mxy}).
\item \eqn{Li} is based on the equations from Li et al. 1993 (\link{Li}).
\item The estimator \eqn{r_{xy}} based on Queller and Goodnight 1989 adapted to pairwise comparisons as described in Oliehoek et al. 2006 (\link{rxy}).
\item \eqn{l_{xy}} is calculated based on Lynch and Ritland 1999 (\link{lxy}, \link{lxy.w}).
\item The estimator \eqn{loiselle} is based on Loiselle et al. 1995 (\link{loiselle}, \link{loiselle.w}).
\item The estimator \eqn{wang.fin} is based on Wang 2002 for a finite sample (\link{wang.fin}, \link{wang.fin.w}).
\item The estimator \eqn{wang} is based on Wang 2002 including bias correction for sample size (\link{wang}, \link{wang.w}).
\item The estimator \eqn{ritland} is based on Ritland 1996 (\link{ritland}, \link{ritland.w}).
\item The estimator \eqn{morans.fin} is based on Hardy and Vekemans 1999 ommiting correction for sample size (\link{morans.fin}, \link{morans.fin.w}).
\item The estimator \eqn{morans} is based on Hardy and Vekemans 1999 with correction for sample size bias (\link{morans}, \link{morans.w}).
}

Detailed information can be found in the references or the descriptions of either subfunction.


For each estimator populations of randomized offsprings and randomized non related individuals are created from a reference population. By default, all empirical data will be used as a reference; however, functions can be forced to use specific data as reference. Be careful when using several populations or specific references. All populations need to be free of any pronounced population structure. References must contain every allele of the tested populations! Any violation of this may lead to strange results or interruptions during calculations.
Based on these reference populations all allele sharing indices and relatedness estimators can be analysed using a logistic regression model. Thresholds for individuals being full-siblings, half-siblings or non-related for each population are calculated. \eqn{\chi^2} statistics are used to calculate whether populations contain more siblings than expected in randomly drawn populations of the same size as the empirical population. Additionally, for different sets of references (full siblings, half siblings and non-related) in comparison to all empirical populations pairwise T-tests are used to calculate whether populations are significantly different in mean relatedness. These results may be biased if the statistical conditions for performing a T-test are not given in the data. You should take care of that before drawing conclusions from these results. Therefore pairwise relatedness for empirical populations as well as randomized populations are passed to the local directory and can be used for downstream analyses. Several plots from each analysis can be exported as .pdf files for visualization of relatedness within populations. For every population geographical locations can be combined with genetic data to find out whether geographic distance has an effect on relatedness. A mantel correlation (Oksanen et al. 2013) plotted via scatterplot including relatedness thresholds from logistic regression are plotted for visualization. 

Additionally, basic information on datasets can be calculated, for example allele and genotype frequencies, \eqn{F_{is}} values within populations or analysis on relative informativity of loci for relatedness analysis. For a quick overview, while omitting statistics, all allele sharing indices and relatedness estimators are calculated for every pairwise combination in one analysis.


}

\details{

\tabular{ll}{
Package: Demerelate\tab \cr
Version: 0.9-1\tab\cr
Date: 2016-08-31\tab\cr
Depends: R (>= 2.15.0), Formula\tab\cr
Imports: fts, mlogit, sfsmisc, vegan\tab\cr
Suggests: MASS\tab\cr
License: GPL (>= 2)\tab\cr
URL: http://www.r-project.org\tab\cr
}

\tabular{ll}{
Demerelate \tab -- Head function to calculate pairwise relatedness\cr
Loci.test \tab -- Function to analyse datasets of diploid genetic information\cr
Emp.calc \tab -- A quick overview of relatedness by omitting statistics\cr
F.stat \tab -- Function to calculate \eqn{F_{is}} values for single populations\cr
demerelpop \tab -- Example data set from randomized diploid markers for three populations\cr
demerelref \tab -- Example data set from randomized diploid markers as reference\cr
demereldist \tab -- Example data set from randomized relative geographic information\cr
}
}

\references{
Armstrong, W. (2012) fts: R interface to tslib (a time series library in c++). by R package version
0.7.7.\cr
Blouin, M., Parsons, M., Lacaille, V. and Lotz, S. (1996) Use of microsatellite loci to classify indi-
viduals by relatedness. Molecular Ecology, 5, 393-401.\cr
Croissant, Y. 2011 mlogit: multinomial logit model. R package version 0.2-2.\cr
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
Maechler, M. and many others. (2012) Utilities from Seminar fuer Statistik ETH Zurich. R package
version 1.0-20.\cr
Nei, M. (1977) F-statistics and analysis of gene diversity in subdivided populations. Annals of Hu-
man Genetics, 41, 225-233.\cr
Nei, M. and Chesser R.K. (1983) Estimation of fixation indices and gene diversities. Annals of
Human Genetics, 47, 253-259.\cr
Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with
a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
Oksanen, J. et al. (2013) vegan: Community Ecology Package. R package version 2.0-8.\cr
Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution,
43, 258-275.\cr
Ritland, K. (1999) Estimators for pairwise relatedness and individual inbreeding coefficients. Ge-
netics Research, 67, 175-185.\cr
Wang, J. (2002) An estimator for pairwise relatedness using molecular markers. Genetics, 160,
1203-1215.\cr
Weir, B.S. and Cockerham, C.C. (1984) Estimating F-Statistics for the analysis of population struc-
ture. Evolution, 38, 1358-1370.\cr
 }

\author{
Philipp Kraemer and Gabriele Gerlach

Maintainer: Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

\keyword{package}
