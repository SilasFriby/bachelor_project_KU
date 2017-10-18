
rm(list=ls()) # clears all

#### OMXC20

data_omxc20 <- read.table('data/omxc20.csv', header = T, sep = ';', dec = ',')
data_omxc20$Date <- as.Date(data_omxc20$Date, format = '%d-%m-%Y')
dim(data_omxc20)
head(data_omxc20)
str(data_omxc20)

plot(1:dim(data_omxc20)[1], data_omxc20$Closing.price)


## log returns

log_return <- log(data_omxc20$Closing.price[-1] / data_omxc20$Closing.price[-dim(data_omxc20)[1]])
hist(log_return)

## sub data

data_omxc20_sub <- subset(data_omxc20, Date > as.Date('2017-01-17'))
dim(data_omxc20_sub)

log_return <- log(data_omxc20_sub$Closing.price[-1] / data_omxc20_sub$Closing.price[-dim(data_omxc20_sub)[1]])
hist(log_return)

## estimate Q-parameters in GBM

dt <- 1/(dim(data_omxc20_sub)[1]/ 5) # one over the number of observations per year
sigma_mle <- sqrt(var(log_return) / dt)
r_mle <- mean(log_return) / dt + sigma_mle / 2

hist(log_return, freq = FALSE)
lines(sort(log_return), dnorm(sort(log_return), mean = (r_mle - sigma_mle / 2) * dt, sd = sigma_mle * sqrt(dt)), lwd = 2)

## test independence of log returns - seems OK

library(stats)

acf(log_return)
Box.test(log_return, lag = 1, type = "Ljung-Box")

## test normal distribution - not OK!

qqnorm(log_return);qqline(log_return, col = 2)


shapiro.test(log_return) 


## option price data


