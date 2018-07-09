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
library(brms)

sc <- make_stancode(bf(score ~ tau + (tau|id/test)),
              data = p, family = gaussian(),
              set_prior('normal(0, 0.1)', class = "b", coef = 'tau')))
