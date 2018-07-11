library(ggplot2)

set.seed(1)
## Data generation model

tau <- 13  # number of time points
J <- 8  # number of tests
alpha_j <- rep(rnorm(J, 0, 1/2), tau)  # subject level performance for test j
gamma_j <- rep(abs(rnorm(J, 0, 1/4)), tau) # subject level sd for test j
bp <- floor(rnorm(1, tau/2, tau/5))
mu_tau <- c(rep(0, J*bp), rep(1, (J*tau - J*bp))) # subject level absence(0)/presence(1) of underlying pathological state
u_0t <- rep(rgamma(tau, 1.5, 10), each = J) # correlated variance (how the subject feels on the day)


X <- data.frame(time = rep(1:tau, each = J),
                test = rep(1:J, tau),
                alpha_j,
                gamma_j,
                mu_tau,
                u_0t,
                bp)

X$score <- -0.1*X$mu_tau*(X$time-X$bp) - X$u_0t + rnorm(J*tau, X$alpha_j, X$gamma_j)

ggplot(X, aes(x = time, y = score, colour = factor(test))) +
  geom_line() +
  geom_point() + 
  geom_vline(aes(xintercept = X$bp[1]), linetype = "dotted")
