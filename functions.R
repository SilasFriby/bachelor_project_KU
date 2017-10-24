

#' Simulate GBM
#'
#' @param x0 - initial value of the process
#' @param n - number of simulation paths
#' @param dt - time step length
#' @param end_time - end time
#' @param mu - drift
#' @param sigma - volatility
#' @param drift - choose to have constant or stochastic drift
#'
#' @return matrix with column representing a sample path

simulate_gmb <- function(x0,n, dt, end_time, mu, sigma, drift = 'constant') {
  
  time_vector <- seq(0, end_time, by = dt)
  m  <- length(time_vector) 
  
  if (drift == 'constant') {
    
    x <- matrix(NA, m, n)
    x[1,] <- x0
    
    for (k in 1:n) {
      for (i in 2:m) {
        dx <- mu * x[i - 1, k] * dt + x[i - 1, k] * sigma * sqrt(dt) * rnorm(1,0,1)
        x[i, k] <- x[i - 1, k] + dx
      }
    }
    
  }
  
  if (drift == 'stochastic') {
    
    x <- matrix(NA, m, n)
    x[1,] <- x0
    
    for (k in 1:n) {
      for (i in 2:m) {
        dx <- mu[i - 1, k] * x[i - 1, k] * dt + x[i - 1, k] * sigma * sqrt(dt) * rnorm(1,0,1)
        x[i, k] <- x[i - 1, k] + dx
      }
    }
    
  }
  
  
  
  
  
  
  
  return(x)
  
}


#' Simulate Vasicek
#'
#' @param x0 - initial value of the process
#' @param n - number of simulation paths
#' @param dt - time step length
#' @param end_time - end time
#' @param kappa - speed of mean reversion
#' @param theta - mean reversion
#' @param sigma - volatility
#'
#' @return matrix with column representing a sample path

simulate_vasicek <- function(x0, n, dt, end_time, kappa, theta, sigma) {
  
  time_vector <- seq(0, end_time, by = dt)
  m  <- length(time_vector) 
  
  x <- matrix(NA, m, n)
  x[1,] <- x0
  
  for (k in 1:n) {
    for (i in 2:m) {
      dx <- kappa * (theta - x[i - 1, k]) * dt + sigma * sqrt(dt) * rnorm(1,0,1)
      x[i, k] <- x[i - 1, k] + dx
    }
  }
  
  return(x)
  
}


vasicek_zcb_price <- 
  function(r0, k, theta, beta, T){
    b.vas <- (1/k)*(1-exp(-T*k)) 
    a.vas <- (theta-beta^2/(2*k^2))*(T-b.vas)+(beta^2)/(4*k)*b.vas^2
    return(exp(-a.vas-b.vas*r0))
  }


VASICEKyield<-function(r,tau,Pparam,riskpremium=0)
{ b<-Pparam[1]+riskpremium
a<-Pparam[2]
sig<-Pparam[3]
Btau<-(1-exp(-a*tau))/a
Atau<-((Btau-tau)*(a^2*b-0.5*sig^2)/a^2 - sig^2*Btau^2/(4*a))
return(r*Btau/tau-Atau/tau)
}


black_scholes_formula <- function(ttm, s, K, r, sigma, dividend) {
  
  d1 <- 1 / (sigma * sqrt(ttm)) * (log(s / K) + (r - dividend + sigma ^ 2 / 2) * ttm)
  d2 <- d1 - sigma * sqrt(ttm)
  
  call_price <- pnorm(d1) * s * exp(-dividend * ttm) - pnorm(d2) * K * exp(-r * ttm) 
  
  return(call_price)
  
}
