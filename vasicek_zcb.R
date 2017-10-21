
VasicekZCBprice <- 
  function(r0, k, theta, beta, T){
    b.vas <- (1/k)*(1-exp(-T*k)) 
    a.vas <- (theta-beta^2/(2*k^2))*(T-b.vas)+(beta^2)/(4*k)*b.vas^2
    return(exp(-a.vas-b.vas*r0))
  }