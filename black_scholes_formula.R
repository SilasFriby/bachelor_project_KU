

black_scholes_formula <- function(t, maturity, s, K, r, sigma) {
  
  d1 <- 1 / (sigma * sqrt(maturity - t)) * (log(s / K) + (r + sigma ^ 2 / 2) * (maturity - t))
  d2 <- d1 - sigma * sqrt(maturity - t)
  
  call_price <- pnorm(d1) * s - pnorm(d2) * K * exp(-r * (maturity - t)) 
  
  return(call_price)
  
}
