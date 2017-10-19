

black_scholes_formula <- function(ttm, s, K, r, sigma, dividend) {
  
  d1 <- 1 / (sigma * sqrt(ttm)) * (log(s / K) + (r - dividend + sigma ^ 2 / 2) * ttm)
  d2 <- d1 - sigma * sqrt(ttm)
  
  call_price <- pnorm(d1) * s * exp(-dividend * ttm) - pnorm(d2) * K * exp(-r * ttm) 
  
  return(call_price)
  
}
