data {
     int<lower=0,upper=1> pre_working_age; // 0 = no, 1 = yes
     int<lower=0,upper=1> working_age; // 0 = no, 1 = yes
     int<lower=0,upper=1> retirement; // 0 = no, 1 = yes
     int<lower=0,upper=1> post_retirement; // 0 = no, 1 = yes
     int<lower=0> N;
     vector[N] x;
     vector[N] y;
}
parameters {
     real<lower=0> alpha1[1*pre_working_age];
     real<lower=0> alpha2[1*working_age];
     real<lower=0> alpha3[1*retirement];
     real<lower=0, upper=1> a1[1*pre_working_age];
     real<lower=0, upper=1> a2[1*working_age];
     real<lower=0, upper=1> a3[1*retirement];
     real<lower=0, upper=1> a4[1*post_retirement];
     real<lower=0> mu2[1*working_age];
     real<lower=0, upper=max(x)> mu3[1*retirement];
     real<lower=0> lambda2[1*working_age];
     real<lower=0> lambda3[1*retirement];
     real<upper=0.05> lambda4[1*post_retirement];
     real<lower=0, upper=1> c;
     real<lower=0> sigma;
}
transformed parameters {
     vector[N] mu_rc;
     vector[N] mu_rc_1;
     vector[N] mu_rc_2;
     vector[N] mu_rc_3;
     vector[N] mu_rc_4;
     vector[N] zero;
     
     for(i in 1:N){
        zero[i] = 0;
     }
     

     mu_rc_1 = pre_working_age==1?a1[1]*exp(-alpha1[1]*x):zero;
     mu_rc_2 = working_age==1?a2[1]*exp(-alpha2[1]*(x - mu2[1]) - exp(-lambda2[1]*(x - mu2[1]))):zero;
     mu_rc_3 = retirement==1?a3[1]*exp(-alpha3[1]*(x - mu3[1]) - exp(-lambda3[1]*(x - mu3[1]))):zero;
     mu_rc_4 = post_retirement==1?a4[1]*exp(lambda4[1]*(x)):zero;
     mu_rc = mu_rc_1 + mu_rc_2 + mu_rc_3 + mu_rc_4 + c;
}
model {
     // likelihood
     y ~ normal(mu_rc, sigma);
     
     //priors
     
     if(pre_working_age==1){
      alpha1 ~ normal(0,1);
      a1 ~ normal(0,0.1);
     }
     if(working_age==1){
      alpha2 ~ normal(0,1);
      a2 ~ normal(0,0.1);
      mu2 ~ normal(25,1);
      lambda2 ~ normal(0,1);
     }
     if(retirement==1){
      alpha3 ~ normal(0,1);
      a3 ~ normal(0,0.1);
      mu3 ~ normal(65,1);
      lambda3 ~ normal(0,1);
     }
     if(post_retirement==1){
      a4 ~ normal(0,0.05);
      lambda4 ~ normal(0,0.01);
     }
     c ~ normal(min(y),0.1);
     sigma ~ normal(0,1);
}
