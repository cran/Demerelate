## Philipp Kraemer function created as part of Demerelate 28.07.2015
# Updated on 06.10.2015

wang.compose <- function(Ps, as)
  
{
  
  b <- .subset2(as,1)
  c <- .subset2(as,2)
  d <- .subset2(as,3)
  e <- .subset2(as,4)
  f <- .subset2(as,5)
  g <- .subset2(as,6)

    # eq. 11
    V <- (1-b)^2*(e^2*f + d*g^2)-(1-b)*(e*f-d*g)^2+2*c*d*f*(1-b)*(g+e)+c^2*d*f*(d+f)
    
    
    # eq. 9
    twogen <- (d*f*((e+g)*(1-b)+c*(d+f))*(Ps[1]-1)+
                 d*(1-b)*(g*(1-b-d)+f*(c+e))*Ps[3]+
                 f*(1-b)*(e*(1-b-f)+d*(c+g))*Ps[2])/V
    
    # eq. 10
    fourgen <- (c*d*f*(e+g)*(Ps[1]+1-2*b)
                + ((1-b)*(f*e^2+d*g^2)-(e*f-d*g)^2)*(Ps[1]-b)
                + c*(d*g-e*f)*(d*Ps[3]-f*Ps[2])-c^2*d*f
                * (Ps[3] + Ps[2] - d - f) - c*(1-b)*(d*g*Ps[3] + e*f*Ps[2]))/V
    
    # equation (1)
    r.basic <- twogen/2 + fourgen
    return(r.basic)
  
}