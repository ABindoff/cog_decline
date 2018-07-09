####  Piecewise regression with informative priors
####  estimate break-point as a model parameter 
####  (rather than as a model selection procedure)
####  incorporate covariates
####  6th July 2018
####  see https://github.com/sambrilleman/2017-Epidemiology
library(reshape2)
library(viridis)
library(purrr)

set.seed(1)


### simulate data, j tests of t time points, one participant "A"
tc <- function(x,t){
    c(rep(0, x[3]), ((x[3]+1)-x[3]):(t-x[3]))*x[4]+ rnorm(t, x[1], sqrt(x[2]))
}

sim <- function(id, j = 9, t = 12, j.mu = 0, j.var = 1){
  s <- data.frame(j.mu = j.mu,
                  j.var = j.var,
                  cp.t = 7,
                  cp.slope = -abs(rnorm(j, 0, 0.5)),
                  covariate = sample(c(1:6), j, replace = TRUE))
  p <- melt(data.frame(apply(s, 1, function(x) tc(x, t))))
  p$tau <- as.double(rep(c(1:t), j))
  p$id <- id
  p$cp <- s$cp.t
  p$covariate <- as.vector(sapply(s$covariate, function(x) rep(x, t)))
  names(p) <- c("test", "score", "tau", "id", "cpT", "covariate")
  p$test <- as.integer(p$test)
  p
}  

j = 9
t = 12
j.mu = rnorm(j, 0, 1)
j.var = rnorm(j, 0, 0.5)^2

  p <- map_dfr(c(1:6), function(x) sim(x, j = j, t = t, j.mu = j.mu, j.var = j.var))


#  fit a piece-wise regression using the known bp and plot the data with the regression slopes
p$cp0 <- ifelse(p$tau > p$cp, p$tau - p$cp, 0)
p$test <- as.integer(p$test)
m0 <- lm(score ~ tau + cp0, p)
p0 <- p
p0$fit <- predict(m0, p)

library(ggplot2)
ggplot(p0, aes(x = tau, y = score, colour = factor(test))) +
  geom_line(size = 3/4) +
  scale_colour_viridis(discrete = TRUE) +
  facet_wrap(~id) +
  geom_line(aes(y = fit), colour = "red", size = 1) +
  geom_vline(data = p0, aes(xintercept = cpT), linetype = "dashed") +
  xlab("Testing phase") +
  scale_x_continuous(breaks = c(1,3,5,7,9,11)) +
  theme(panel.background = element_blank(),
        axis.line = element_line()) 

# fire up rstan
library(rstan)

stancode.filepath <- "stancode_randomchangecorr_randomtestcorr.stan"

standata <- list(N = length(p$score),
                 Npat = length(unique(p$id)),
                 Ntest = length(unique(p$test)),
                 y = p$score,
                 age = p$tau,
                 id = p$id,
                 test = p$test,
                 betakp_lower = min(p$tau),
                 betakp_upper = max(p$tau),
                 Npred = length(p$score),
                 Npat_pred = 2,
                 #Ntest_pred = 9,
                 age_pred = p$tau,
                 id_pred = p$id,
                 #test_pred = p$test,
                 zeros4 = c(0,0,0,0))
  

stanmonitor <- c("beta", 
                 "betakp", 
                 "y_sd", 
                 "u_1_sd", 
                 "u_1_Sigma", 
                 "u_1_Corr",
                 "u_2_sd",
                 "u_2_Corr",
                 "y_mu_pred", 
                 "y_pred", 
                 "alpha_tosave",
                 #"alpha_tosave2",
                 "lp__")

stanfit <- stan(file    = stancode.filepath, 
                data    = standata, 
                pars    = stanmonitor, 
                chains  = 1, 
                cores   = 3, 
                iter    = 3000, 
                warmup  = 1000,
                control = list(adapt_delta = 0.85, max_treedepth = 15))

#=======================================
# Some diagnostics for the fitted model
#=======================================

# Trace plots of the MCMC samples for various model parameters
traceplot(stanfit, 
          pars = c("beta", "betakp"), 
          inc_warmup = TRUE)
traceplot(stanfit, 
          pars = c("u_sd"), 
          inc_warmup = TRUE)
traceplot(stanfit, 
          pars = c("u_Corr"), 
          inc_warmup = TRUE)
traceplot(stanfit, 
          pars = c("lp__"), 
          inc_warmup = TRUE)

# Diagnostics for Hamiltonian Monte Carlo sampler 
# Notes:
#   See page 15 of the RStan vignette (http://mc-stan.org/interfaces/rstan)
#   for some discussion of these diagnostics
summary(do.call(rbind, 
                args = get_sampler_params(stanfit, inc_warmup = FALSE)),
        digits = 2)


#=============================
# Summary of model parameters
#=============================

print(stanfit, pars = c("beta", 
                        "betakp", 
                        "y_sd", 
                        "u_1_sd", 
                        "u_1_Sigma", 
                        "u_1_Corr"))

#=====================================================================
# Plot observed height measurements with fitted trajectories overlaid
#=====================================================================

# Extract MCMC samples for the model
mcmcsamples.pred <- extract(stanfit)

# Calculate the mean trajectory based on the posterior mean for each of the model parameters
pred.data <- data.frame(cbind(id = standata$id_pred, 
                              age = standata$age_pred))
pred.alphas <- data.frame(cbind(id = c(1:length(unique(standata$id_pred))), 
                                apply(mcmcsamples.pred$alpha_tosave, 
                                      MARGIN = c(2,3), 
                                      FUN = mean)))
pred.data <- merge(pred.data, pred.alphas, by = "id")
pred.data[, "y.predtraj"] <- pred.data$V2 + 
  pred.data$V3 * pmin(pred.data$age - pred.data$V5, 0) + 
  pred.data$V4 * pmax(pred.data$age - pred.data$V5, 0)

# Calculate percentiles of the posterior predictive distribution
pred.p025 <- data.frame(cbind(id = standata$id_pred, 
                              age = standata$age_pred, 
                              y.predlimit = apply(mcmcsamples.pred$y_pred, 
                                                  MARGIN = 2, 
                                                  FUN = quantile, probs = .025) ))
pred.p975 <- data.frame(cbind(id = standata$id_pred, 
                              age = standata$age_pred, 
                              y.predlimit = apply(mcmcsamples.pred$y_pred, 
                                                  MARGIN = 2, 
                                                  FUN = quantile, probs = .975) ))
pred.limits <- rbind(pred.p975, pred.p025[order(pred.p025$age, decreasing = TRUE), ])

# Observed data
  obs.data <- data.frame(cbind(id = standata$id, 
                             age = standata$age, 
                             y = standata$y))
  library(MASS)    # package for mvrnorm
  library(doBy)    # package for summaryBy
  library(rstan)   # package for running Stan from within R
  library(lattice) # package for xyplot
  library(latticeExtra)
# Create plot
a <- xyplot(y.predlimit ~ age | id, 
            data = pred.limits[pred.limits$id %in% c(1:10), ],
            panel = panel.polygon,
            xlab = "Test phase", ylab = "z-score", strip=FALSE, 
            scales = list(alternating = FALSE, relation = "sliced"),
            layout = c(5,2), col = 'lightgray', border = 'lightgray')
b <- xyplot(y ~ age | id, 
            data = obs.data[obs.data$id %in% c(1:10), ],
            panel = panel.xyplot,
            type = 'p', pch = 19, col = 'black')
c <- xyplot(y.predtraj ~ age | id,
            data = pred.data[pred.data$id %in% c(1:10), ],
            panel = panel.xyplot,
            type = 'l', lty = 2, col = 'black')
print(a + as.layer(b) + as.layer(c))
