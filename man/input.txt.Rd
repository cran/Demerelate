% File src/library/base/man/input.txt.Rd
     \name{input.txt}
     \alias{input.txt}
     \title{Reads different types of tables and returns an object.}
     \description{
       Internal function to prepare inputdata for Demerelate.
     }
     \usage{
     input.txt(tab.txt, mod)
     }
     \arguments{
        \item{tab.txt}{Can be either a filename which will be passed to an object or an object in the R workspace}
        \item{mod}{Can be either "dist", "pop" or "ref.pop". Different informations are printed}
     }

      \value{
		  \item{tab}{Object containing dataframe as inputdata.}
		  
	      }


\author{Philipp Kraemer, <philipp.kraemer@uni-oldenburg.de>}

     \seealso{
       \code{
	     \link{Demerelate}
	     \link{inputformat}
	      }
     }
     \examples{
      ## internal function not intended for direct use
    
     
     }
