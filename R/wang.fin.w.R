### P. Kraemer 13.5.2015
## Basic functions from Wang 2002
# Updated 06.10.2015 - new name wang.w.fin

# wang.fin calculates from finite populations ommiting formula 12-14 to estimate from a sample
## Basic functions from Wang 2002

wang.fin.w <- function(allele.column, ref.pop = NA)
{
  

p <- as.numeric(as.character(unlist(ref.pop[allele.column])))

# p=vector of observed gene frequencies i,j,k .. n 
# m=given by category of P

    a2 <- sum(p^2)
    
    # eq. 13  
    a3 <- sum(p^3)
    
    # eq. 14
    a4 <- sum(p^4)
    
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
