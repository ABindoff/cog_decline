####  Piecewise regression with informative priors
####  estimate break-point as a model parameter 
####  (rather than as a model selection procedure)
####  4th July 2018
library(reshape2)
library(viridis)

set.seed(1)


### simulate data, j tests of t time points, one participant "A"
tc <- function(x,t){
  c(rep(0, x[3]), ((x[3]+1)-x[3]):(t-x[3]))*x[4]+ rnorm(t, x[1], sqrt(x[2]))
}

sim <- function(id, j = 9, t = 12, j.mu = 0, j.var = 1){
  s <- data.frame(j.mu = j.mu,
                  j.var = j.var,
                  cp.t = 7,
                  cp.slope = -abs(rnorm(j, 0, 0.66)))
  p <- melt(data.frame(apply(s, 1, function(x) tc(x, t))))
  p$tau <- as.double(rep(c(1:t), j))
  p$id <- id
  p$cp <- s$cp.t
  names(p) <- c("test", "score", "tau", "id", "cpT")
  p
}  

j = 9
t = 12
j.mu = rnorm(j, 0, 1)
j.var = rnorm(j, 0, 0.5)^2

p <- sim("A", j = j, t = t, j.mu = j.mu, j.var = j.var)


#  fit a piece-wise regression using the known bp and plot the data with the regression slopes
p$cp0 <- ifelse(p$tau > p$cp, p$tau - p$cp, 0)
m0 <- lm(score ~ tau + cp0, p)
p0 <- p
p0$fit <- predict(m0, p)

library(ggplot2)
ggplot(p0, aes(x = tau, y = score, colour = test)) +
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

### build a model
model <- '
data { 
int<lower=1> N;  // number of obsvn
real score[N];      
real tau[N];      
}

parameters { 
real intercept;                 // predicted score at breakpoint
real beta_1;              // healthy cognitive state
real beta_2;               // cognitive decline
real<lower = 1, upper = 11> bp; // bp_tau
real<lower = 0> error;          // sd
} 

transformed parameters{
vector[N] conditional_mean; // the estimated average score at each tau
real beta_diff;

beta_diff = beta_2 - beta_1;  

// conditional_mean dependent on whether tau is before or after bp
for (i in 1:N) {
if (tau[i] < bp) {
conditional_mean[i] = intercept + beta_1 * (tau[i] - bp);
} else {
conditional_mean[i] = intercept + beta_2 * (tau[i] - bp);
}
}
}

model {
// Set priors
intercept ~ normal(0, 2);  // Average score at breakpoint
beta_1 ~ normal(0, 1);  // Slope at healthy cognitive function
beta_2 ~ normal(-1, 1);   // Slope after cognitive decline
bp ~ normal(6, 6);           // Breakpoint tau
error ~ normal(0, 3);        // Residual error

// assumed data generation model
for (i in 1:N) {
score[i] ~ normal(conditional_mean[i], error);
}
}

generated quantities {
vector[N] sim_score;               // Simulate new scores using estimated params
vector[N] log_lik;               // probably need this
vector[12] sim_conditional_mean; // pretty pics

// Compute conditional means for tau
for (i in 1:12) {
if (i < bp) {
sim_conditional_mean[i] = intercept + beta_1 * (i - bp);
} else {
sim_conditional_mean[i] = intercept + beta_2 * (i - bp);
}
}

for (i in 1:N) {
sim_score[i] = normal_rng(conditional_mean[i], error);
log_lik[i] = normal_lpdf(score[i] | conditional_mean[i], error);
}
}
'

tau <- p$tau
score <- p$score

data <- list(
  tau = tau,
  score = score,
  N = length(score)
)


###  fit model using stan
m1 <- stan(model_code = model, 
                   data = data)

### check results
print(m1,
      par = c("intercept", "bp", "beta_1", "beta_2", "beta_diff", "error"))

plot(m1, pars = c("bp", "beta_1", "beta_2", "intercept"))
