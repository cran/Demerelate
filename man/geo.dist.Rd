% File src/library/base/man/geo.dist.Rd
     \name{geo.dist}
     \alias{geo.dist}
     \title{Calculates geographic distances.}
     \description{
       Internal function to prepare geographic distances for the combined analysis with genetic relatedness.
     }
     \usage{
     geo.dist(pop1, pop2, onlypairs = FALSE, value)
     }
     \arguments{
        \item{pop1}{Specific type of dataframe as in \code{inputformat}. Population one used for calculations. Individuals passed to rows of resulting matrix. Inputformat should be standart with x and y coordinate mode for Demerelate.}
       \item{pop2}{Specific type of dataframe as in \code{inputformat}. Population two used for calculations. Individuals passed to columns of resulting matrix. Inputformat should be standard with x and y coordinate mode for Demerelate.}
       \item{onlypairs}{If set as TRUE geographic distances is calculated only for diagonal comparisons in matrix. If set as FALSE lower triangel is calculated.}
       \item{value}{String defining method to calculate geographic distances. Can be set as "relative" or "decimal".}
     }

      \details{Two different methods of distance calculations are implemented in Demerelate. If using "relative", distances will be calculated from x-y coordinates using normal pythagoratic mathematics. When working with geographic positions value needs to be set to "decimal". x and y coordinates need to be given as geographic positions in decimal degrees.}

      \value{
		  \item{matrix.share}{Object containing geographic distances.}
		  
	      }

\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{Demerelate}
	     \link{allele.sharing}
	     \link{inputformat}
	      }
     }
     \examples{
     ## internal function not intended for direct use
    
     
     }
