### P. Kraemer 13.5.2015
# Updated 06.10.2015
## Basic functions from Wang 2002

wang.w <- function(allele.column, ref.pop = NA)
{
  
N <- as.numeric(names(ref.pop[allele.column]))
p <- ref.pop[[allele.column]]


# p=vector of observed gene frequencies i,j,k .. n 
# N=Sample size
# m=given by category of P 

    a2 <- (N*(sum(p^2))-1)/(N-1)
    
    # eq. 13  
    a3 <- (N^2*sum(p^3)-3*(N-1)*a2-1)/((N-1)*(N-2))
    
    # eq. 14
    a4 <- (N^3*sum(p^4)-6*(N-1)*(N-2)*a3-7*(N-1)*a2-1)/(N^3-6*N^2+11*N-6)
    
    a5 <- a2^2

# Make short variables / (2*as[1]-as[2]) == u

u <- 2*a2-a3

b <- (2*a5-a4)/(u)

c <- (a2-2*a5+a4)/(u)

d <- (4*(a3-a4))/(u)

e <- (2*(a2 - 3*a3 + 2*a4))/(u)

f <- (4*(a2-a5-2*a3+2*a4))/(u)

g <- (1-7*a2+4*a5+10*a3-8*a4)/(u)


return(c(b,c,d,e,f,g,u))

}
