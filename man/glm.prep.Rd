% File src/library/base/man/glm.prep.Rd
     \name{glm.prep}
     \alias{glm.prep}
     \title{Preparation for multinomial logistic regression.}
     \description{
       Function combines randomized full and half siblings with empirical values and non related random individuals to one multinomial logistic regression to calculate relatedness thresholds.
     }
     \usage{
     glm.prep(empirical.list, offfull.list, 
              offhalf.list, offnon.list)
     }
     \arguments{
		  \item{empirical.list}{Mean relatedness for empirical population over all loci. Population information from inputdata are ommited.} 
		  \item{offfull.list}{Mean relatedness of populations randomized as full siblings over all loci from the reference population.}
		  \item{offhalf.list}{Mean relatedness of populations randomized as half siblings over all loci from the reference population.} 
		  \item{offnon.list}{Mean relatedness of populations of randomized non relatives over all loci from the reference population.} 
     }

      \details{The function uses the package mlogit from Croissant 2011 to combine all the information from reference populations to one multilogistic regression model using mlogit(..). Thresholds are calculated and used for downstream calculations.}

      \value{
		  \item{sumlmr}{Summary of mlogit model}
		  \item{half}{Thresholds for full and half siblings}
		  
	      }

\references{
      Croissant, Y. 2011 mlogit: multinomial logit model R package version 0.2-2.}

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{offspring}
	     \link{random.pairs}
	      }
     }
     \examples{
     ## internal function not intended for direct use
    
     
     }
