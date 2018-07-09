library(reshape2)
library(viridis)
library(purrr)
library(rstan)

set.seed(1)


### simulate data, j tests of t time points, one participant "A"
tc <- function(x,t){
  c(rep(0, x[3]), ((x[3]+1)-x[3]):(t-x[3]))*x[4]+ rnorm(t, x[1], sqrt(x[2]))
}

sim <- function(id, j = 9, t = 12, j.mu = 0, j.var = 1){
  s <- data.frame(j.mu = j.mu,
                  j.var = j.var,
                  cp.t = 7,
                  cp.slope = -abs(rnorm(j, 0, 0.33)),
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


# int N; //the number of observations
# int J; //the number of groups
# int K; //number of columns in the model matrix
# int id[N]; //vector of group indeces
# matrix[N,K] X; //the model matrix
# vector[N] y; //the response variable
# 

# fire up rstan
library(rstan)

stancode.filepath <- "simple_hierarchical_piecewise_model.stan"

p0 <- p
p0$x <- p0$tau
p0$y <- p0$score
X <- model.matrix(~ x + y, p0)


# test is random effect, subject is ignored
standata <- list(N = length(p$score),
                 J = length(unique(p$test)),
                 K = ncol(X),  # intercept, slope
                 id = p$test,
                 X = X,
                 y = p$score)
                 
stancode.filepath = "simple_hierarchical_piecewise_model.stan"


stanfit <- stan(file    = stancode.filepath, 
                data    = standata, 
                chains  = 1, 
                cores   = 3, 
                iter    = 3000, 
                warmup  = 500,
                control = list(adapt_delta = 0.85, max_treedepth = 15))

mcmc_hier<-extract(stanfit)


traceplot(stanfit, 
          pars = c("gamma"), 
          inc_warmup = TRUE)



#plot average response to explanatory variables
X_new<-model.matrix(~x+y,data=data.frame(x=seq(-2,2,by=0.2),y=0))
#get predicted values for each MCMC sample
pred_x1<-apply(mcmc_hier$gamma,1,function(beta) X_new %*% beta)
#now get median and 95% credible intervals
pred_x1<-apply(pred_x1,1,quantile,probs=c(0.025,0.5,0.975))
#same stuff for the second explanatory variables
X_new<-model.matrix(~x+y,data=data.frame(x=0,y=seq(-2,2,by=0.2)))
pred_x2<-apply(mcmc_hier$gamma,1,function(beta) X_new %*% beta)
pred_x2<-apply(pred_x2,1,quantile,probs=c(0.025,0.5,0.975))


cols<-viridis(10)
par(mfrow=c(1,2),mar=c(4,4,0,1),oma=c(0,0,3,5))
plot(y~X[,2],pch=16,xlab="Temperature",ylab="Response variable",col=cols[id])
lines(seq(-2,2,by=0.2),pred_x1[1,],lty=2,col="red")
lines(seq(-2,2,by=0.2),pred_x1[2,],lty=1,lwd=3,col="blue")
lines(seq(-2,2,by=0.2),pred_x1[3,],lty=2,col="red")
plot(y~X[,3],pch=16,xlab="Nitrogen concentration",ylab="Response variable",col=cols[id])
lines(seq(-2,2,by=0.2),pred_x2[1,],lty=2,col="red")
lines(seq(-2,2,by=0.2),pred_x2[2,],lty=1,lwd=3,col="blue")
lines(seq(-2,2,by=0.2),pred_x2[3,],lty=2,col="red")
mtext(text = "Population-level response to the two\nexplanatory variables with 95% CrI",side = 3,line = 0,outer=TRUE)
legend(x=2.1,y=10,legend=paste("Gr",1:10),ncol = 1,col=cols,pch=16,bty="n",xpd=NA,title = "Group\nID")
