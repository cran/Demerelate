\name{Similarity.Index}
\alias{Bxy}
\alias{Sxy}
\alias{Mxy}
\alias{rxy}
\alias{loiselle}
\alias{loiselle.w}
\alias{lxy}
\alias{lxy.w}
\alias{Li}
\alias{morans}
\alias{morans.w}
\alias{morans.fin}
\alias{morans.fin.w}
\alias{ritland}
\alias{ritland.w}
\alias{wang}
\alias{wang.w}
\alias{wang.fin}
\alias{wang.fin.w}
\alias{wang.compose}

\title{
Pairwise genetic similarity calculated based on different relatedness estimators or allele shraing indices
}
\description{
Pairwise estimates are based on one of several inidices:

Allele/genotype sharing indices, which count the number of shared alleles in different ways:\cr
\enumerate{
\item \eqn{B_{xy}} (number of alleles shared) as described in Li and Horvitz 1953.
\item \eqn{S_{xy}} (proportion of shared alleles) as described in Lynch 1988.
\item \eqn{M_{xy}} (genotype sharing) as described in Blouin et al. 1996.\cr}

Relatedness estimators, which account for the similarity in allele composition of two individuals by chance (identity by state; IBS) based on reference allele frequencies:
\enumerate{
\item The estimator \eqn{r_{xy}} based on Queller and Goodnight 1989 adapted to pairwise comparisons as described in Oliehoek et al. 2006.
\item \eqn{l_{xy}} is calculated based on Lynch and Ritland 1999.
\item The estimator \eqn{ritland} is based on Ritland 1996.
\item The estimator \eqn{wang.fin} is based on Wang 2002 for a finite sample.
\item The estimator \eqn{morans.fin} is based on Hardy and Vekemans 1999 omitting correction for sample size.\cr}

Relatedness estimators, which additionally unbias for sample size effects:
\enumerate{
\item \eqn{Li} is based on the equations from Li et al. 1993.
\item The estimator \eqn{loiselle} is based on Loiselle et al. 1995.
\item The estimator \eqn{wang} is based on Wang 2002 including bias correction for sample size.
\item The estimator \eqn{morans} is based on Hardy and Vekemans 1999 with correction for sample size bias.\cr
}

Single locus similarities are either simply averaged over loci for each pairwise comparison (\eqn{B_{xy}}, \eqn{S_{xy}}, \eqn{M_{xy}}), weighted for each locus before averaging (\eqn{r_{xy}}, \eqn{l_{xy}}, \eqn{ritland}, \eqn{Li}) or the multilocus estimate is weigthed by the average of weigths for each pairwise comparison over loci (\eqn{wang.fin}, \eqn{morans.fin}, \eqn{loiselle}, \eqn{wang}, \eqn{morans}).
}

\usage{

Bxy(row, data, pop1, pop2, allele.column, ref.pop = NA)
Sxy(row, data, pop1, pop2, allele.column, ref.pop = NA)
Mxy(row, data, pop1, pop2, allele.column, ref.pop = NA)
rxy(row, data, pop1, pop2, allele.column, ref.pop = NA)
Li(row, data, pop1, pop2, allele.column, ref.pop = NA)
ritland(row, data, pop1, pop2, allele.column, ref.pop = NA)
ritland.w(row, data, pop1, pop2, allele.column, ref.pop = NA)

lxy(row, data, pop1, pop2, allele.column, ref.pop = NA)
lxy.w(row, data, pop1, pop2, allele.column, ref.pop = NA)

loiselle(row, data, pop1, pop2, allele.column, ref.pop = NA)
loiselle.w(row, data, pop1, pop2, allele.column, ref.pop = NA)

wang.fin(row, data, pop1, pop2, allele.column, ref.pop = NA)
wang.fin.w(allele.column, ref.pop = NA)
wang(row, data, pop1, pop2, allele.column, ref.pop = NA)
wang.w(allele.column, ref.pop = NA)
wang.compose(Ps, as)

morans.fin(row, data, pop1, pop2, allele.column, ref.pop = NA)
morans.fin.w(row, data, pop1, pop2, allele.column, ref.pop = NA)

morans(row, data, pop1, pop2, allele.column, ref.pop = NA)
morans.w(row, data, pop1, pop2, allele.column, ref.pop = NA)


}

\arguments{
  \item{row}{A numeric value which sets the row of \code{data} used for calculations}
  \item{data}{A dataframe with all pairwise combinations of pop1 and pop2, column one set the name, column two the row number of the individual in pop1 and column three set the row number of individual two in pop2.}
  \item{pop1}{Specific dataframe of type \link{inputformat}. Population one used for calculations. Inputformat needs to be standard three digits mode for Demerelate.}
  \item{pop2}{Specific dataframe of type \link{inputformat}. Population two used for calculations. Inputformat needs to be standard three digits mode for Demerelate.}
  \item{allele.column}{Numeric value - It equals the number of of the first column in the dataframe containing allele information. Order of loci in both populations needs to be exactly equal.}
  \item{ref.pop}{Reference population used for relatedness calculations.}
  \item{Ps}{Vector of observed similarity classes P as a part of the twogen and fourgen similarity as described in Wang 2002}
  \item{as}{Vector of parameters to unbias the observed similarity class of P based on observed gene frequencies. Parameters are equvalent to those denoted as b-g and u in Wang 2002}
}

\details{
\eqn{B_{xy}}\cr
The similarity is calculated based on the average number of alleles shared between two individuals. If there are at least three alleles in the locus - A, B and C, two diploid individuals may have four different states of similarity. If all alleles are the same in individuals Bxy=1, if two are the same between individuals (either AB vs AB or AA vs AB) Bxy=0.5, if only one allele is the same (AB vs AC) Bxy=0.25 and if no allele is shared Bxy=0 (Li and Horvitz 1953).

\eqn{S_{xy}}\cr
The similarity is calculated based on the average number of allele positions, which share the same allele in both individuals. If there are at least three alleles in the locus - A, B and C, two diploid individuals may have four different states of similarity. If all alleles are the same in individuals - Sxy=1. If both individuals are heterozygous and both alleles are present in both individuals (AB vs AB) Sxy=1. If one individual is homozygous for a shared allele (eg. AA vs AB) Sxy=0.75. If only one allele is the same in both individuals (AB vs AC) Sxy=0.5 and if no allele is shared Sxy=0 (Lynch 1988).

\eqn{M_{xy}}\cr
Sharing rate is calculated according to shared allele positions i.e. 0, 1 or 2 shared allele positions for diploids. A sharing rate of 0 is calculated if no alleles are shared, a rate of 0.5 if only one allele position is equal in individuals (AC vs AB or AA vs AB) and a rate of 1 if individuals match in both allele positions (AA vs AA or AB vs AB) (Blouin et al. 1996).

\eqn{r_{xy}}\cr
The estimator \eqn{r_{xy}} based on Queller and Goodnight 1989 adapted to pairwise comparisons as described in Oliehoek et al. 2006 is calculated as follows:
\tabular{ll}{
 \tab \cr
 \tab \eqn{r_{xy,l} = \frac{0.5(I_{ac}+I_{ad}+I_{bc}+I_{bd})-p_{a}-p_{b}}{1+I_{ab}-p_{a}-p_{b}}}\cr
 \tab \eqn{I_{a-d} = } allele identities of individual I in locus l\cr
 \tab \eqn{p_{a-b} = } frequencies of allele a or b in reference populations\cr
}
The equation is consulted twice for each individual pairing via a RE-RAT procedure. This means individuals are switched by calculating \eqn{r_{xy}} and \eqn{r_{yx}} and averaging for the final pairwise estimate.

\eqn{l_{xy}}\cr
In Lynch and Ritland 1999 \eqn{l_{xy}} is referred to as \eqn{r_{xy}} arcoding to equation (5a). For multilocus estimates weights over loci are calculated according to equation (6a and 7a) (Lynch and Ritland 1999).  
\tabular{ll}{
 \tab \cr
 \tab \eqn{l_{xy,l} = \frac{p_{a}*(S_{bc}+S_{bd})+p_{b}*(S_{ac}+S_{ad})-4p_{a}p_{b}}{(1+S_{ab})*(p_{a}+p_{b})-4p_{a}p_{b}}}\cr
 
 \tab \eqn{w_{xy,l} = \frac{(1+S_{ab})*(p_{a}+p_{b})-4p_{a}p_{b}}{2p_{a}p_{b}}}\cr
 \tab \eqn{S_{a-d} = } allele identities of individual S in locus l\cr
 \tab \eqn{p_{a-b} = } frequencies of allele a or b in reference populations\cr
}
Each pairwise estimate is weighted with \eqn{w_{xy,l}} for each locus and calculated via RE-RAT as described for \code{rxy}. The final pairwise estimate is divided by the sum of weights for loci \eqn{W_{xy}}, which is averaged analogous to RE-RAT procedure of the pairwise relatedness estimate. 

\eqn{Li}\cr
The estimator is calculated according to Li et al. 1993 equation 9 corrected for the average similarity for unrelated individuals based on reference allele frequencies as \eqn{U=2a_{2}+a_{3}}. The inital similarity \eqn{S_{xy}} is calculated based on Lynch 1988.

\tabular{ll}{
 \tab \cr
 \tab \eqn{Li_{xy,l} = \frac{S-U}{1-U}}\cr

\tab \eqn{S = } similarity after Lynch 1988\cr
 \tab \eqn{U = S(unrelated) = 2a_{2}+a_{3} } after Li et al. 1993\cr
}

\eqn{ritland}\cr
Ritlands original estimator from Ritland 1996 (equation 5) is calculated according to Lynch and Ritland 1999 by multiplying the final estimate with 2. As the sum of all p is equal to 1 the equation is simplyfied to:

\tabular{ll}{
 \tab \cr
 \tab \eqn{ritland_{xy,l} = \frac{2}{n-1}[(\sum\frac{S_{i}}{p_{i}})-1]}\cr

\tab \eqn{S_{i} = } allele identity in both individuals\cr
 \tab \eqn{p_{i} = } reference frequency of allele i\cr
  \tab \eqn{n = } number of different alleles in locus l\cr
}

The single locus estimate is averaged over loci. The basic similarity is calculated for allele i as 0.25 if i is present in both individuals, 0.5 if i is present in both and one individual is homozogious for i and 1 if both individuals are homozygous for i. Note that this similarity measure is equivalent to Blouin's approach \eqn{M_{xy}} (Blouin et al. 1996) if summed over all alleles for each pairwise comparison.

\eqn{loiselle}\cr
The estimator first described in Loiselle et al. 1995 is implemented as described by Hardy and
Vekemans 2015. The frequency of each allele in individuals (i.e. 0.5 or 1 for diploids) in a pairwise
comparison are combined and corrected for the allele frequency in the reference population. The
product of corrected shared allele frequencies is additionally corrected for sample size bias and
combined over loci via weighting with the polymorphic index (\eqn{\sum{p_{i}} * (1 - p_{i})}).

\tabular{ll}{
 \tab \cr
 \tab \eqn{loiselle_{xy,l} = \sum{
				\frac{
					\sum{(p_{ila}-p{la})*(p_{jla}-p_{la})}
				       +\sum{(p_{la}*(1-p_{la}))}
                                      }
                                {n_{l}-1}}/
			     \sum{\sum{p_{la}*(1-p_{la})}}
           }\cr

 \tab \eqn{p_{ila} = } frequency of allele a in individual i of locus l
\cr
  \tab \eqn{n = } number of different alleles in locus l\cr
}

\eqn{wang}\cr
Equations described by Wang 2002 are implemented as follows: 

A binary vector \eqn{P} of classes of observed similarities is calculated for each pairwise comparison. One of the 4 categories is set as 1, all remaining P are set as 0. \eqn{P_{1}} is 1 if both animals are homozygous for the same allele or both are heterozygous with the same allele combination. \eqn{P_{2}} is 1 if one individual is homozygous and the other is heterozygous sharing only 1 allele. \eqn{P_{3}} is 1 if only 1 allele is shared between individuals with 1 copy per individual. All other comibnations fall into the category four \eqn{P_{4}}. Note that the categories follow the general classification as described for \eqn{S_{xy}} described by Lynch 1988.

The estimator descirbed by Wang 2002 is relatively complex and can be discussed only superficial here. The probability of each category i.e. the joint probability of genotypes can be estimated using the two-gene (\eqn{\Theta}) and four-gene (\eqn{\Delta}) coefficient based on the sum of powers of allele frequencies (\eqn{a_{m}=\sum{p_{i}^{m}}}) summing up each probability of each category of \eqn{P}. Formulas used for calculations can be found in Wang 2002 Equation (9), (10) and (10). Combining these estimates into r=\eqn{\frac{\Theta}{2}+\Delta} yields the estimate by Wang (2002). By using the estimator \code{wang} each expected sum of powers of allele frequency is corrected for sample size N (\eqn{\bar{a}_{2}}, \eqn{\bar{a}_{3}}, \eqn{\bar{a}_{4}}) compare eqation (12), (13) and (14) in Wang 2002. 

Finally, to get an multilocus estimate each single locus estimate is corrected for the average similarity value for unrelated individuals \eqn{u=2a_{2}+a_{3}} as described in Li et al. 1993. Each single term of the estimate \eqn{\Theta} and \eqn{\Delta}, namely \eqn{b-g} and each \eqn{P} are unbiased by \eqn{u} and the combined estimate is weigthed with the sum of \eqn{u} over loci.


\eqn{wang.fin}\cr
This estimate is based on the same equations as \code{wang}. However, for the use of finite samples the bias correction of sample size (N) is omitted. Instead of \eqn{\bar{a}_{2}}, \eqn{\bar{a}_{3}}, \eqn{\bar{a}_{4}} ((12), (13), and (14) in Wang 2002) the pure sum of powers of allele frequencies (\eqn{a_{m}=\sum{p_{i}^{m}}}) are used.

\eqn{morans}\cr
The estimator morans refers to morans I, which is widely used as estimator for spatial autocorrela-
tion. It is described in Hardy and Vekemans 1999 as estimator for genetic relatedness. The approach
is similar to loiselle (Loiselle et al. 1995) or \eqn{r_{xy}} (Queller and Goodnight 1989) by correcting each
individual allele frequency by the allele frequency of the reference sample for each shared allele of
each pairwise comparison. Each shared allele is then weighted by the variance of individual allele
frequencies and a term for sample size bias. In order to calculate the estimator over loci the sum of
pairwise estimates over loci is weighted by the sum of weights over loci for each comparison.


\tabular{ll}{
 \tab \cr
 \tab \eqn{morans_{xy,l} = \frac{
				\sum{
					(p_{ila}-p_{la})*(p_{jla}-p_{la})
				    }
                                }
                            {\sum{Var(p_{ila}+1/(n_{l}-1))}}
           }\cr

\tab \eqn{p_{ila} = } frequency of allele a in individual i of locus l
\cr
 \tab \eqn{n = } number of different alleles in locus l
\cr
}

\eqn{morans.fin}\cr
The estimator morans.fin refers to the same calculations as morans but omitting weights for sample
size bias. Thus it should only be applied for finite samples (Hardy and Vekemans 1999).


}

\value{
The value of similarity is returned according to the arguments defining the individuals compared.
}

\references{
Armstrong, W. (2012) fts: R interface to tslib (a time series library in c++). by R package version
0.7.7.\cr
Blouin, M., Parsons, M., Lacaille, V. and Lotz, S. (1996) Use of microsatellite loci to classify indi-
viduals by relatedness. Molecular Ecology, 5, 393-401.\cr
oissant, Y. 2011 mlogit: multinomial logit model. R package version 0.2-2.\cr
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
ture. Evolution, 38, 1358-1370.
}

\author{
Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>
}

\seealso{
\code{\link{allele.sharing}
      \link{Demerelate}
      \link{Emp.calc}
	}
}

\examples{
## internal function not intended for direct use
}
