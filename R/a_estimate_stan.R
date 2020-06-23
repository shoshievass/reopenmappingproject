
# install.packages("rstan")
# install.packages("gridExtra")
# install.packages("outbreaks")

rm(list = ls())

## set up directory
dir<-Sys.getenv("HOME")

if (dir=="/Users/hanyang") {
  proj <- paste(dir, "Documents","GitHub","reopenmappingproject", sep="/")
  stopifnot(dir.exists(proj))
  dataPath <- paste(proj, "data", sep="/")
  codePath <- paste(proj, "R",    sep="/")
  outPath  <- paste(proj, "output", sep="/")
  parmPath <- paste(proj, "parameter", sep="/")
} 
#check if project directory is properly set up
stopifnot(endsWith(proj, "reopenmappingproject"))
setwd(codePath)



library(rstan)
library(gridExtra)
library(outbreaks)
rstan_options (auto_write = TRUE)
# options (mc.cores = parallel::detectCores ())
options (mc.cores = 1)
set.seed(3) # for reproductibility

## set up global varibales and functions
source(paste(codePath,"2_programs.R",sep='/'))
source(getCodePath("1_config.R"))

# load contact matrix from synthetic population

contact <- "_W1-S1-N1-E2-M2"
m  <- "1600"
Cmat<<-loadData(place, contact)





### load covid death and cases data
COVID <-checkLoad(deathData)
covid<-COVID[COVID$msa==m,]

# total population
N <-mean(covid$population)

# times series of cases
cases <- covid$cases

# times
n_days <- length(cases)
t <- seq(0, n_days, by = 1)
t0 <- 0 
t <- t[-1]

#initial conditions
i0 <- 5
s0 <- N - i0
y0 = c(S = s0, E=0, I = i0, R = 0)


# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases, epsilon=0.3)

# number of MCMC steps
niter <- 2000

# compile model
model <- stan_model("sir_negbin.stan")

# estimate
fit_sir_negbin <- sampling(model,
                           data = data_sir,
                           iter = niter,
                           chains = 4)

# show estimation results
pars=c('beta', 'gamma', "R0", "recovery_time")
print(fit_sir_negbin, pars = pars)

stan_dens(fit_sir_negbin, pars = pars, separate_chains = TRUE)


# plot predicted number of infected
smr_pred <- cbind(as.data.frame(summary(
  fit_sir_negbin, pars = "pred_cases", probs = c(0.05, 0.5, 0.95))$summary), t, cases)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.)) + 
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Number of infected")

