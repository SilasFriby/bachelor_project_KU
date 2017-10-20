
## denne kode er "stjålet" fra Rolf

rm(list = ls())

#### useful functions

VASICEKnegloglike<-function(param,data,times)
{ n.obs<-length(data)
dum<-sort.list(times)
data<-data[dum]
times<-times[dum]
delta<-diff(times,1)
mv<-data[1:(n.obs-1)]*exp(-param[2]*delta)+param[1]*(1-exp(-param[2]*delta))
variance<-param[3]^2*(1-exp(-2*param[2]*delta))/(2*param[2])
VASICEKnegloglike<--sum(log(dnorm(data[2:n.obs],mv,sqrt(variance))))/n.obs

}

VASICEKyield<-function(r,tau,Pparam,riskpremium=0)
{ b<-Pparam[1]+riskpremium
a<-Pparam[2]
sig<-Pparam[3]
Btau<-(1-exp(-a*tau))/a
Atau<-((Btau-tau)*(a^2*b-0.5*sig^2)/a^2 - sig^2*Btau^2/(4*a))
VASICEKyield<-r*Btau/tau-Atau/tau
}



#### Read the data

data_us_rates <- fread("data/us_rates_2001_2017.csv")
data_us_rates <- as.data.frame(data_us_rates)
data_us_rates[, -1] <- data_us_rates[, -1] / 100 # rates in %
str(data_us_rates)
data_us_rates <- data_us_rates[complete.cases(data_us_rates), ] # remove NA's
head(data_us_rates)

# subset data
data_us_rates <- subset(data_us_rates, date >= as.Date('2010-01-01'))

obs<-data_us_rates[,2]
N<-length(obs)-1
dt<-1/250
data<-obs[2:(N+1)]
lagdata<-obs[1:N]



#### estimate P-parameters from observations using maximum likelihood

## Closed form estimators

bhat<-(sum(data*lagdata) - sum(data)*sum(lagdata)/N)/(sum(lagdata*lagdata) - sum(lagdata)*sum(lagdata)/N)

kappahat<--log(bhat)/dt

ahat<-sum(data)/N-bhat*sum(lagdata)/N
thetahat<-ahat/(1-bhat)
s2hat<-sum((data-lagdata*bhat-ahat)^2)/N

sigmahat<-sqrt(2*kappahat*s2hat/(1-bhat^2))

## plot the data (or more accurately, the short rate)

plot(obs,type='l', xlab="date", ylab="US 3M rate", main="US short rate 1952-2004")
abline(h=thetahat, lty=3)

## parameters found by numerical optimization - just to see if they match the closed-form estimate

# mle_vasicek<-optim(par=c(thetahat,kappahat,sigmahat),fn=VASICEKnegloglike,method = "BFGS",  data=obs,times=dates)

mle_vasicek <- c(thetahat,kappahat,sigmahat)

#### Estimate Q mean reversion parameter from observed yield curves 

# idea: yield curves are used in pricing, hence the the vasicek yield curve (with Q-parameters) should match the observed yield curves


## analyze average yield curves

maturities <- c(1/12, 1/4, 1/2, 1, 2, 3, 5, 7, 10, 20, 30)
meanyield<-2:dim(data_us_rates)[2]

for(i in 1:(dim(data_us_rates)[2] - 1)) meanyield[i]<-mean(data_us_rates[,i+1])

# plot average observed yield curve
plot(maturities,meanyield,type='l', ylab="mean/average/typical yield", main="Typical yield curves w/ US  2001-2017 data")


# if P=Q, ie. riskpremium = lamda= 0, these curves should 
# match the average yields well

points(maturities,VASICEKyield(meanyield[1],maturities,mle_vasicek,riskpremium=0),type='l',lty=2)
# text(maturities[6],VASICEKyield(meanyield[1],maturities[6],mle_vasicek,riskpremium=0)+0.001, "Typical Vasicek model curve, I (r0 = avr. short rate, risk premium=0) ")

# clearly, they don't. So let's try with a simple extension: A constant 
# riskpremuim, dW^Q = dW^P + lambda dt, or differently 
# theta^Q = theta^P + tilde{lambda} 
#          = (theta^P  +what the code calls the riskpremium) 

# let's estimate tilde{lambda} as what gives the best fit of the 
# average (over time) shape of the yield curve

rpfit<-function(rp){
  rpfit<-sum(abs(meanyield-VASICEKyield(meanyield[1],maturities,mle_vasicek,riskpremium=rp)))
}

rp<-optim(par=c(0),fn=rpfit,method = "BFGS")$par

theory<-VASICEKyield(meanyield[1],maturities,mle_vasicek,riskpremium=rp)
points(maturities,theory,type='l',lty=2)
# text(maturities[6],theory[6], "Typical Vasicek model curve, III (r0=avr. short rate, risk premium='what fits best') ")

# theory says that the model should match the yield curve not 
# just "on average" (that would be a "necessary condition")
# but EVERY day

plot(maturities,data_us_rates[(N+1),2:12],type='l',ylab="yield", main="On the last day in the sample")
# text(maturities[6],data_us_rates[(N+1),7], "Observed yield curve")


points(maturities,VASICEKyield(data_us_rates[(N+1),2],maturities,mle_vasicek,riskpremium=rp),type='l',lty=2)
# text(maturities[6],VASICEKyield(data_us_rates[(N+1),2],maturities[6],mle_vasicek,riskpremium=rp), "Vasicek model 'prediction' (w/ all the estimates) ")

# AND that's why we calibrate & 'back to Bjork'
# Silas: we calibrate in order for the model to match observe yield today (because a trusted model should at least fit today's reality right??)