// adapted from https://github.com/sambrilleman/2017-Epidemiology
data {
  int<lower=0> N;            // number of observations
  int<lower=1> Npat;         // number of individuals  
  int<lower=1> Ntest;         // number of tests
  real y[N];        // outcome data
  real<lower=1> age[N];      // explanatory variable (age) data
  int<lower=1,upper=Npat> id[N]; // patient id for each observation
  int<lower=1,upper=Ntest> test[N];  //test id for each observation
  vector<lower=0,upper=0>[4] zeros4; // mean vector for random effects distribution
  int<lower=1> betakp_lower; // lower bound for (prior on) mean knot point
  int<lower=1> betakp_upper; // upper bound for (prior on) mean knot point
  
  int<lower=0> Npred;        // number of predicted observations
  int<lower=0> Npat_pred;    // number of patients to predict observations for
  real<lower=0> age_pred[Npred]; // explanatory variable (age) for prediction
  //int<lower=0> Ntest_pred;         // number of tests to predict observations for 
  int<lower=1,upper=Npat> id_pred[Npred]; // patient id for each predicted observation
  //int<lower=1,upper=Ntest> test_pred[Ntest]; // test id for each predicted observation
}

parameters {
  vector[3] beta;            // fixed effects, intercept and slopes
  real<lower=betakp_lower,upper=betakp_upper> betakp;
                             // fixed effects, knotpoint (bounding specified to help convergence)
  vector<lower=0>[4] u_1_sd;   // level 2 error sd (sds of the random effects u_1[j])
  vector<lower=0>[4] u_2_sd;  // level 3 error sd  (sd of random effect test u_2[j])
  real<lower=0> y_sd;        // level 1 error sd
  vector[4] u_1[Npat];         // random effects (level 2) errors
  vector[4] u_2[Ntest];       // re test errors
  cholesky_factor_corr[4] L_u_1_Corr;
                             // cholesky factor for the random effects correlation matrix
  cholesky_factor_corr[4] L_u_2_Corr;// cholesky factor for random effects corr mat tests 
}

transformed parameters {  
  vector[4] alpha[Npat];     // random effects
  vector[4] alpha_2[Ntest];  // random effects tests
  real y_mu[N];              // mean parameter based on regression equation

  //==========================
  // calculate random effects
  //==========================
  
  for (i in 1:Npat) {
    for (k in 1:3) alpha[i,k] = beta[k] + u_1[i,k];
    alpha[i,4] = betakp + u_1[i,4];
  }
  
  for (i in 1:Ntest) {
    for (k in 1:3) alpha_2[i,k] = beta[k] + u_2[i,k];
    alpha_2[i,4] = betakp + u_2[i,4];
  }
  
  //=====================
  // regression equation
  //=====================
  
  for (j in 1:N) {
    if (age[j] < alpha[id[j],4]) 
      y_mu[j] = alpha[id[j],1] + alpha[id[j],2] * (age[j] - alpha[id[j],4]) + alpha_2[id[j],1] + alpha_2[id[j],2] * (age[j] - alpha_2[id[j],4]);
    else  
      y_mu[j] = alpha[id[j],1] + alpha[id[j],3] * (age[j] - alpha[id[j],4]) + alpha_2[id[j],1] + alpha_2[id[j],3] * (age[j] - alpha_2[id[j],4]);     
  } 
}

model {

  //========
  // priors
  //========
  
  beta[1] ~ normal(0, 4); // prior: fixed effect, intercept
  beta[2] ~ normal(0, 1);   // prior: fixed effect, slope before knot
  beta[3] ~ normal(-1, 1);   // prior: fixed effect, slope after knot
  betakp ~ uniform(betakp_lower,betakp_upper);
                            // prior: fixed effect, knot point

  u_1_sd[1] ~ cauchy(0,5);    // prior: random effect sd, intercept
  u_1_sd[2] ~ cauchy(0,5);    // prior: random effect sd, slope before knot
  u_1_sd[3] ~ cauchy(0,5);    // prior: random effect sd, slope after knot
  u_1_sd[4] ~ cauchy(0,5);    // prior: random effect sd, knot point
  
  u_2_sd[1] ~ cauchy(0,5);
  u_2_sd[2] ~ cauchy(0,5);
  u_2_sd[3] ~ cauchy(0,5);
  u_2_sd[4] ~ cauchy(0,5);
  
  
  y_sd ~ cauchy(0,5);       // prior: level 1 error sd
  
  L_u_1_Corr ~ lkj_corr_cholesky(1);
  L_u_2_Corr ~ lkj_corr_cholesky(1);

               // prior: cholesky factor for random effects correlation matrix
               // NB. this prior is the "lkj correlation distribution" with shape parameter 1 
               // which is equivalent to a uniform distribution over the possible correlation 
               // matrices (where a shape parameter > 1 would have resulted in an upside down
               // U-shaped distribution with the mode being located at the identity matrix)
  //=============================
  // random effects distribution
  //=============================
  
  for (i in 1:Npat) u_1[i] ~ multi_normal_cholesky(zeros4, diag_pre_multiply(u_1_sd, L_u_1_Corr));
  for (i in 1:Ntest) u_2[i] ~ multi_normal_cholesky(zeros4, diag_pre_multiply(u_2_sd, L_u_2_Corr));

                             // NB. the second parameter here is the cholesky factor L 
                             // (for the correlation matrix). It only uses the sd rather 
                             // than the variances since Sigma = L*L'
  
  //==================
  // model likelihood
  //==================
  
  y ~ normal(y_mu, y_sd); // likelihood for the observed data
  
}

generated quantities {
  real y_pred[Npred];      // predicted outcome
  real y_mu_pred[Npred];   // predicted mean
  corr_matrix[4] u_1_Corr;   // random effects correlation matrix
  matrix[4,4] u_1_Sigma;     // random effects covariance matrix
  vector[4] alpha_tosave[Npat_pred];
                           // monitor random effects for a subset of patients only
                           // (for plotting predictions) and do not monitor 'alpha' 
                           // in the model above (since it consumes too much memory!)
  corr_matrix[4] u_2_Corr;   // random effects correlation matrix
  matrix[4,4] u_2_Sigma;     // random effects covariance matrix
  //vector[4] alpha_tosave_2[Ntest_pred];
  
  //==================================================
  // predicted mean outcome using regression equation
  //==================================================

  for (i in 1:Npat_pred) {  
    alpha_tosave[i] = alpha[i];
  }
  
  // for (i in 1:Ntest_pred) {  
  //   alpha_tosave_2[i] = alpha_2[i];
  // }
  
  for (j in 1:Npred) {
    if (age_pred[j] < alpha[id_pred[j],4]) 
      y_mu_pred[j] = alpha[id_pred[j],1] + alpha[id_pred[j],2] * (age_pred[j] - alpha[id_pred[j],4]);
    else  
      y_mu_pred[j] = alpha[id_pred[j],1] + alpha[id_pred[j],3] * (age_pred[j] - alpha[id_pred[j],4]);      
  
    y_pred[j] = normal_rng(y_mu_pred[j], y_sd);
  }
  

  
  //=====================================================
  // recover the correlation and covariance matrices
  // using the cholesky factor of the correlation matrix
  //=====================================================
  
  u_1_Corr = multiply_lower_tri_self_transpose(L_u_1_Corr);    
	            // correlation matrix: u_Corr = L_u_1_Corr * L_u_1_Corr'
  
	u_1_Sigma = quad_form_diag(u_1_Corr, u_1_sd);
	            // covariance matrix: u_Sigma = diag(u_1_sd) * u_Corr * diag(u_1_sd)
  u_2_Corr = multiply_lower_tri_self_transpose(L_u_2_Corr);    
	            // correlation matrix: u_Corr = L_u_2_Corr * L_u_2_Corr'
  
	u_2_Sigma = quad_form_diag(u_2_Corr, u_2_sd);
	            // covariance matrix: u_Sigma = diag(u_2_sd) * u_Corr * diag(u_2_sd) 
}









