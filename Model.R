library(ggplot2)
library(rstan)

#install.packages("rstan", type = "source")
setwd("/Users/xiaofanliang/Documents/Minerva/Academic/CS146/Final")

dataset <- read.csv("weekly_in_situ_co2_mlo.csv")
colnames(dataset) <- c("date", "CO2")

#Change the date formate 
dates <- as.Date(dataset$date,
                 format = "%m/%d/%y")

form_date <- as.Date(ifelse(dates > Sys.Date(), format(dates, "19%y-%m-%d"), format(dates)))

dataset$date <- form_date

#Plot and show the original data that need to be modelled. 
plot(dataset$date[1:250],dataset$CO2[1:250], type='l')
#lines(dataset$date[1:250], dataset$CO2[1:250], type='l')
#plot(dataset$date,dataset$CO2)

#----------------------------------------------------------------------------------------#
###########TRY GENERATE FUTURE DATA AND CONFIDENCE INTERVAL FOR FUTURE DATA#########
###########FAILED. GENERATED DATA IS A POSTERIOR DISTRIBUTION. HOW CAN I VISUALIZE THAT?###

#Try it with the model professor provided 

f <- function(dt, c0, c1, c2, c3){
  return(c0 + c1 * dt  + c2 * cos(2 * pi * dt / 365.25  + c3))
}

seqnum <- seq(0,3042,length=3042)

n <- 250
n_future <- 250
ppm <- dataset$CO2[1:n]

#date <- seqnum

data <- list(
  n = n,
  n_future = n_future,
  ppm = ppm[1:n] #actual data from 1:n 
)

#int date;

stan_model <- "
data {
  int n;        //number of datapoints 
  int n_future; //number of future datapoints generated 
  real ppm[n];  //CO2 ppm data, each datapoint should be a real value
}

parameters {
  real c0;          //intercept 
  real c1;          //trend coefficient 
  real c2_prime; 
  real phase_x;
  real phase_y;
  real sigma_prime;
}

transformed parameters {
  real<lower=0> c2;
  real<lower=0,upper=2*pi()> phase;
  real<lower=0> sigma;

  c2 = exp(c2_prime);
  phase = atan2(phase_x, phase_y)/ (2*pi());
  sigma = exp(sigma_prime);
}

model {
  c0 ~ normal(315,15);
  c1 ~ normal(0.03, 2);
  c2_prime ~ normal(0, 10);
  phase_x ~ normal(0, 1);
  phase_y ~ normal(0,1);
  sigma_prime ~ normal(5, 5);

  for(t in 1:n) {
  ppm[t] ~ normal(c0 + c1 * t + c2 * cos(2 * pi() * t / 365.25 + phase), sigma);
  }
}

generated quantities {
  real x_future[n_future];

  for(t in 1:n_future) {
  x_future[t] = normal_rng(c0 + c1 * (t+n) + c2 * cos(2 * pi() * (t+n) / 365.25 + phase), sigma);
  }
}
"

#where the hack does x_future come in?? 

#Fit the Stan model. This will take about 2 minutes.
fit <- stan(
  model_code = stan_model,
  data = data,
  # number of Markov chains
  chains = 4,             
  # number of warmup iterations per chain
  warmup = 1000,          
  # total number of iterations per chain
  iter = 2000,            
  # number of cores (using 2 just for the vignette)
  cores = 1,              
  # show progress every 'refresh' iterations
  refresh = 1000,         
  control = list(adapt_delta = 0.999)
)

print(fit, par=c('c0', 'c1', 'c2', 'phase','sigma'), probs=c(.05, .5, 0.95))
samples <- extract(fit)

# Plot 
results <- apply(samples$x_future, 2, quantile, probs=c(0.025, 0.975))  # Compute quantiles of the predicted values

plot(
  1:n, ppm,
  col='black', type='l',
  xlim=c(1, n+n_future),
  ylim=c(min(c(results[1,], ppm)), max(c(results[2,], ppm))),
  main='Data, Future Data, and predicted 95% interval')

lines((n+1):(n+n_future), results[1,], lty='dashed', col='blue')   # 95% interval in blue
lines((n+1):(n+n_future), results[2,], lty='dashed', col='blue')
abline(v=n, col='red')

######---------Other ways of examing how good the samples are---------################ 
acf(samples$c0, main="Autocorrelation of c0 samples")
acf(samples$c1, main="Autocorrelation of c1 samples")
acf(samples$c2, main="Autocorrelation of c2 samples")
acf(samples$phase, main="Autocorrelation of c3 samples")
acf(samples$sigma, main="Autocorrelation of sigma samples")

# Plot histograms of the parameter samples
pairs(fit, pars = c("c0", "c1", "c2", "phase","sigma"))

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
###### Try giving a confidence interval predicted to the existing data #####################
setwd("/Users/xiaofanliang/Documents/Minerva/Academic/CS146/Final")

dataset <- read.csv("co2_mlo.csv")
colnames(dataset) <- c("date", "CO2", "Timestep")

rescale <- function(x, min, max){
  # Rescales data vector: 0 < data < 1
  ((x - min) / (max - min))
}

inverse_rescale <- function(x, original.max, original.min) {
  # Will transform ppm / predictions back to original scale
  x * (original.max - original.min) + original.min
}

#date <- seqnum
#//what parameter means. What is a reasonable range according to data
#//c2 = [0,20] -> c2’
#//c2’ = log(c2); pick a dist so cdf at 95% is [-2.3, 3]
#//N(0.35, 1.3)-> -2.3 to be 2std away from mean.
#//Draw samples from prior to see check what the prior means

n <- 250
ppm <- rescale(dataset$CO2[1:n], min(dataset$CO2[1:n]), max(dataset$CO2[1:n]))
t <- dataset$Timestep[1:n]
n_future <- 250

data <- list(
  n = n,
  ppm = ppm, 
  t = t,
  n_future = n_future,
  t_future = t[n:n_future]
)

#int date;

stan_model <- "
data {
  int n;            //number of datapoints 
  int n_predict;    //number of predicted datapoints generated (based on model provided pro.)
  real t[n];        //Timestep transformed from dates 
  real ppm[n];      //CO2 ppm data, each datapoint should be a real value
  int n_future;     //number of future timestep
  real t_future[n_future]; // future timesteps (test and real future)

parameters {
  real c0_prime;    //intercept 
  real c1_prime;    //trend coefficient 
  real c2_prime;    //amplitude 
  real phase_prime; //phase shift 
  real sigma_prime; //variance 
}

transformed parameters {
  real<lower=0> c0;
  real<lower=0> c1;
  real<lower=0> c2;
  real<lower=0> phase;
  real<lower=0> sigma;

  c0 = exp(c0_prime);
  c1 = exp(c1_prime);
  c2 = exp(c2_prime);
  phase = exp(phase_prime);
  sigma = exp(sigma_prime);
}

model {
  //Priors 
  c0_prime ~ normal(0,1);
  c1_prime ~ normal(0,1);
  c2_prime ~ normal(0,1);
  phase_prime ~ normal(0,1);
  sigma_prime ~ normal(0,1);

  //Likelihood 
  for(n in 1:n) {
    ppm[n] ~ normal(c0 + c1 * t[n] + c2 * cos(2 * pi() * t[n] / 365.25 + phase), sigma);
    }
}

generated quantities {
  real ppm_future[n_future];

  for(n in 1:n_predict) {
    ppm_future[n] = normal_rng(c0 + c1 * t_future[n] + c2 * cos(2 * pi() * t_future[n] / 365.25 + phase), sigma);
  }
}
"

#where the hack does x_future come in?? 

#Fit the Stan model. This will take about 2 minutes.
fit <- stan(
  model_code = stan_model,
  data = data,
  # number of Markov chains
  chains = 4,             
  # number of warmup iterations per chain
  warmup = 1000,          
  # total number of iterations per chain
  iter = 2000,            
  # number of cores (using 2 just for the vignette)
  cores = 1,              
  # show progress every 'refresh' iterations
  refresh = 1000,         
  control = list(adapt_delta = 0.999)
)

print(fit, par=c('c0', 'c1', 'c2', 'phase','sigma'), probs=c(.05, .5, 0.95))
samples <- extract(fit)

# Plot 
results <- apply(samples$x_predict, 2, quantile, probs=c(0.025, 0.975))  # Compute quantiles of the predicted values

plot(
  1:n, ppm,
  col='black', type='l',
  xlim=c(1, n),
  ylim=c(min(c(results[1,], ppm)), max(c(results[2,], ppm))),
  main='Data, Future Data, and predicted 95% interval')

lines(1:n, results[1,], lty='dashed', col='blue')   # 95% interval in blue
lines(1:n, results[2,], lty='dashed', col='blue')

abline(v=n, col='red')

#samples$x_predict