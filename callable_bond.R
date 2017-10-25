rm(list = ls())

source('functions.R')
source('present_value_functions.R')
library(reshape2)
library(ggplot2)


rho_temp <- 1
r_range <- seq(0,10)
pv <- rep(NA, length(r_range))


for (r_temp in r_range) {
  
r_index <- which(r_range == r_temp)

## initialize

# simulation
n  <- 10 ^ 2 # monte carlo simulation paths 
end_time  <- 30
dt_r <- 1/240
time_vector_r <- seq(0, end_time, by = dt_r)

dt <- 1/4
time_vector <- seq(0, end_time, by = dt)
m <- length(time_vector)

# short rate parameters
r0 <- r_temp/100 #0.04
kappa <- 0.01
theta_Q <- 0.05
sigma_r <- 0.01


## simulate short rate under Q

short_rate <- simulate_vasicek(x0 = r0,
                               n = n,
                               dt = dt_r,
                               end_time = end_time,
                               kappa = kappa,
                               theta = theta_Q,
                               sigma = sigma_r)

# ## plot short rate paths
# 
# data <- data.frame('time' = time_vector_r, short_rate[, 1])
# data <- melt(data,  id = c('time'))
# 
# ggplot(data, aes(time, value)) +
#   geom_line(aes(colour = variable))


#### MC scheme for callable bond

v <- matrix(NA, m, n)

# prob of PP for exogeneuos reasons
lambda <- 0
p_e <- 1 - exp(-lambda * dt)

# prob of PP for both rate and exogeneuos reasons
rho <- rho_temp
p_r <- 1 - exp(-(rho + lambda) * dt)

# coupon rate annual
q <- 0.01

# face value
face_value_initial <- 100
face_value <- function(t, q, face_value_initial, C) {
  (1 + q) ^ t * face_value_initial - C * sum((1 + q) ^ (0:(t - 1)))
}

# coupon
C <- face_value_initial * q / (1 - (1 + q) ^ -m)

# transaction costs (as a percentage of the remaining principal)
X <- 0 / 100


## boundary condition at maturity

boundary <- C
v[m, ] <- boundary
test <- matrix(NA, m - 1, n)

## test optimality condition

for (i in m:2) {
  
  int <- colSums(short_rate[which(time_vector_r == time_vector[i - 1]):which(time_vector_r == time_vector[i]), ]) * dt_r
  discount_factor <- exp(-int)
  temp_face_value <- face_value(i - 1, q, face_value_initial, C)
  test[i - 1, ] <- discount_factor * v[i,] > temp_face_value * (1 + X)
  v[i - 1, ] <- ifelse(test[i - 1, ],
                       (1 - p_r) * discount_factor * v[i,] + p_r * temp_face_value + C,
                       (1 - p_e) * discount_factor * v[i,] + p_e * temp_face_value + C)
  
  
}

# ## loop backwards trough time 
# 
# test <- matrix(NA, m - 1, n)
# for (i in m:2) {
#   
#   int <- sum(short_rate[which(time_vector_r == time_vector[i - 1]):which(time_vector_r == time_vector[i]), ]) * dt_r
#   temp_face_value <- face_value(i - 1, q, face_value_initial, C)
#   discount_factor <- exp(-int)
#   
#   if (i == m) { # no pp on the last step
#     
#     boundary <- discount_factor * rep(C, n) # boundary
#     v[i - 1, ] <- boundary + C  
#     
#   } else {
#     
#     test[i - 1, ] <- discount_factor * v[i, ] > temp_face_value * (1 + X) # optimality condition
#     v[i - 1, ] <- ifelse(test[i - 1, ],
#                          temp_face_value + C, # X not include seen from the investor's perspective
#                          discount_factor * v[i, ] + C) # X not include seen from the investor's perspective
#     
#     
#   }
#   
# }


# ## plot callable bond paths
# 
# data <- data.frame('time' = time_vector, v)
# data <- melt(data,  id = c('time'))
# 
# ggplot(data, aes(time, value)) +
#   geom_line(aes(colour = variable))
# 
# hist(v[1,])

pv[r_index] <- mean(v[1, ])

}

pv
plot(r_range, pv, type = 'l')
abline(h = face_value_initial * (1 + q))
