data {
     int<lower=0> N;
     vector[N] x;
     vector[N] y;
     }
     parameters {
     real<lower=0> alpha1;
     real<lower=0> alpha2;
     real<lower=0, upper=1> a1;
     real<lower=0, upper=1> a2;
     real<lower=0, upper=1> a4;
     real<lower=0, upper=1> c;
     real<lower=0> mu2;
     real<lower=0> lambda2;
     real<lower=0> lambda4;
     real<lower=0> sigma;
     }
     
     transformed parameters {
     vector[N] mu_rc;
     
     mu_rc = a1*exp(-alpha1*x) + a2*exp(-alpha2*(x - mu2) - exp(-lambda2*(x - mu2))) + a4*exp(lambda4*(x)) + c;
     }
     
     model {
     // likelihood
     y ~ normal(mu_rc, sigma);
     
     //priors
     alpha1 ~ normal(0,1);
     alpha2 ~ normal(0,1);
     a1 ~ normal(0,1);
     a2 ~ normal(0,1);
     a4 ~ normal(0,1);
     c ~ normal(0,1);
     mu2 ~ normal(25,1);
     lambda2 ~ normal(0,1);
     lambda4 ~ normal(0,0.1);
     sigma ~ normal(0,1);
     }
