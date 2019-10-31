setwd("~/OneDrive/R/ESG/publish")
library(readxl)  # allows reading excel files; assumes there is a header row
library("loo")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# last two lines are some recommended Stan options

ts <- (c(1,2,3,5,7,10,20,30)) #the maturities
opy = 1461/28   # observations per year, 1461 weeks in 28 years here
y <- as.matrix(read_excel('us_01_05_18_to_10_04_19.xlsx')) 
#rates by week, with header row
N <- nrow(y)
U <- ncol(y)
c(N,U)  #will feed this to Stan so Stan code can be more general:
df = list(N=N,U=U,y=y,ts=ts,opy=opy)

set.seed(6)  #Stan with some options:
us_1 <- stan(file = 'CIR.stan', data = df, verbose = FALSE, chains = 4, iter = 2000, 
             control = list(adapt_delta = 0.86, max_treedepth = 17))

set.seed(6)  #same for ode 
us_1 <- stan(file = 'CIR_ode.stan', data = df, verbose = FALSE, chains = 4, iter = 2000, 
             control = list(adapt_delta = 0.86, max_treedepth = 17))

#compute loo after Stan run
log_spread_1 <- extract_log_lik(us_1) 
loo_spread_1 <- loo(log_spread_1)
loo_spread_1
#can get loo for each chain 
loo(log_spread_1[1:500,])
loo(log_spread_1[501:1000,])
loo(log_spread_1[1001:1500,])
loo(log_spread_1[1501:2000,])

#write out means by chain
out_us_1 <- get_posterior_mean(us_1)
write.csv(out_us_1, file="out_CIR.csv")

#useful to print and plot some parameter info
paras <- c("kaprn", "omrn", "beta", "kap")
plot(us_1, pars = paras, show_density = TRUE, ci_level = 0.8, fill_color = "purple")
plot(us_1, plotfun = "hist", pars = paras)
print(us_1, pars=paras, probs=c(.025, 0.2, 0.5, 0.8, 0.975), digits_summary = 5)

#extract the samples to have parameter distributions
us_1_ss = extract(us_1, permuted = FALSE) # this gets all the samples
#Need permuted = FALSE to get it in array form
dim(us_1_ss) # shows dimensions, like 500 x 4 x 846; samples, chains , elements
#elements includes everything calculated, including likelihoods;
#usually will want to discard some as file gets big
#can collapse samples and chains to get two dimensions
dim(us1_ss) = c(2000,846)
#but if more than 1 chain, doing this loses variable names
#tricky, as matrices come out in different order than in posterior mean file
