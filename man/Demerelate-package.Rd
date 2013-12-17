\name{Demerelate-package}
\alias{Demerelate-package}
\docType{package}
\title{Demerelate --- Algorithms to estimate pairwise relatedness within populations based on allele sharing}

\description{
The package Demerelate provides several functions to calculate relatedness of individuals based on diploid genetic markers. Following indices can be calculated:

\enumerate{
\item \eqn{B_{xy}} (number of alleles shared) as described in Li and Horvitz 1953.
\item  \eqn{M_{xy}} (genotype sharing) as described in Blouin et al. 1996. Sharing rate is calculated according to shared allele positions i.e. 0, 1 or 2 shared allele postions for diploids. A sharing rate of 0 is calculated if no alleles are shared, a rate of 0.5 if one allele position is equal in individuals and a rate of 1 if individuals match in both allele positions.
\item The estimator \eqn{r_{xy}} based on Queller and Goodnight 1989 adapted to pairwise comparisons as described in Oliehoek et al. 2006 is calculated as follows:
\tabular{ll}{
 \tab \cr
 \tab \eqn{r_{xy,l} = \frac{0.5(I_{ac}+I_{ad}+I_{bc}+I_{bd})-p_{a}-p_{b}}{1+I_{ab}-p_{a}-p_{b}}}\cr
 \tab \eqn{I_{a-d} =} allele identities of individual I in locus l\cr
 \tab \eqn{p_{a-b} =} frequencies of allele a or b in reference populations\cr}
}

For each estimator randomized offsprings and randomized non related populations are created from a reference population. By default as reference all populations will be used, however functions can be forced to use specific data as reference. Be careful when using several populations or specific references. All populations need to be free of any pronounced population structure. References must contain every allele of the tested populations! Any violation of this may lead to strange results. Additionally, geographic distances can be used combined with genetic data to find out whether geographic distance has an effect on relatedness.

Allele sharing indices as \eqn{B_{xy}} or \eqn{M_{xy}} are analysed using a logistic regression model. Thresholds for individuals being full-siblings, half-siblings or non-related for each population are calculated. \eqn{\chi^2} statistics are used to calculate whether populations contain more siblings than expected in randomly drawn populations of the same number of individuals as the empirical population. The pairwise estimator \eqn{r_{xy}} is calculated for different sets of references (full siblings, half siblings and non-related). T-tests are used to calculate whether significantly more siblings are present in populations than expected. Relatedness matrices for empirical populations as well as randomized populations are passed to the local directory and can be used for downstream analyses. Several plots from each analysis are exported as .pdf files for visualization of relatedness within populations. Additionally, basic information on datasets can be calculated for example allele and genotype frequencies, \eqn{F_{is}} values within populations or analysis on relative informativity of loci for relatedness analysis.



}

\details{

\tabular{ll}{
Package: Demerelate\tab \cr
Version: 0.8-1\tab\cr
Date: 2013-12-17\tab\cr
Depends: R (>= 2.15.0), Formula\tab\cr
Imports: fts, mlogit, sfsmisc, vegan\tab\cr
Suggests: MASS\tab\cr
License: GPL (>= 2)\tab\cr
URL: http://www.r-project.org\tab\cr
}

\tabular{ll}{
Demerelate \tab -- Head function to calculate pairwise relatedness\cr
Loci.test \tab -- Function to analyse datasets of diploid genetic information\cr
F.stat \tab -- Function to calculate \eqn{F_{is}} values for single populations\cr
demerelpop \tab -- Example data set from randomized diploid markers for three populations\cr
demerelref \tab -- Example data set from randomized diploid markers as reference\cr
demereldist \tab -- Example data set from randomized relative geographic information\cr
}
}

\references{
Armstrong, W. (2012) fts: R interface to tslib (a time series library in c++). by R package version 0.7.7.\cr
Blouin, M.S. et al. (1996) Use of microsatellite loci to classify individuals by relatedness. Molecular Ecology, 5, 393-401.\cr
Croissant, Y. 2011 mlogit: multinomial logit model. R package version 0.2-2.\cr
Li C.C. and Horvitz D.G. (1953) Some methods of estimating the inbreeding coefficient. American Journal of Human Genetics 5, 107-17.\cr
Maechler, M. and myn others. (2012) Utilities from Seminar fuer Statistik ETH Zurich. R package version 1.0-20.\cr
Nei, M. (1977) F-statistics and analysis of gene diversity in subdivided populations. Annals of Human Genetics, 41, 225-233.\cr
Nei, M. and Chesser R.K. (1983) Estimation of fixation indices and gene diversities. Annals of Human Genetics, 47, 253-259.\cr
Oliehoek, P. A. et al. (2006) Estimating relatedness between individuals in general populations with a focus on their use in conservation programs. Genetics, 173, 483-496.\cr
Oksanen, J. et al. (2013) vegan: Community Ecology Package. R package version 2.0-8.\cr
Queller, D.C. and Goodnight, K.F. (1989) Estimating relatedness using genetic markers. Evolution, 43, 258-275.\cr
Weir, B.S. and Cockerham, C.C. (1984) Estimating F-Statistics for the analysis of population structure. Evolution, 38, 1358-1370.
 }

\author{
Philipp Kraemer and Gabriele Gerlach

Maintainer: Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

\keyword{package}
