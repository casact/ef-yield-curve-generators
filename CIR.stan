// Stan code for a single CIR process
data { //these must be open in R and are read in
  int<lower=0> N;    //# rows 
  int<lower=0> U;    //# columns 
  matrix[N,U] y;     // the data to fit 
  real ts[U];       // taus = maturities
  real opy;          //observations per year
}
transformed data {//to output the log likelihoods to do  loo
  int z[U*N,2];  //index for rows and columns of each point
  int o;        //# observations in y
  o = 0;
  for (n in 1:N) {
    for (t in 1:U) {
        o = o+1;
        z[o,1] = n;
        z[o,2] = t;
    }
  }
}

parameters {  //unless specified in model section, priors are uniform on the ranges
 real<lower=0, upper=1> kaprn; //risk-neutral kappa
 real<lower=0, upper=1>  omrn;  //risk-neutral omega   
 real<lower=0, upper=1> beta;    
 real<lower=0, upper=1> kap;   //real world kappa
  real<lower=-6, upper=-1> log_sigma_y; //log stdev of residuals 
  vector<lower=0, upper=20>[N] rC;  //fitted CIR short rate for each period
}  

transformed parameters {//variables used in calculations have to be declared
 real om;
 vector[U] tau;
 matrix[1,U] C;   //C(tau). Matrix with one row for later generalization
 matrix[1,U] D;   //D(tau)
 vector[U] W; //W = 1/Q for CIR
 real sigma_y;      
 real h;  //for CIR C and D
 vector[N-1] Cmean; //all the means and variances of the CIR process after the first
 vector[N-1] Cvar;  //priors are distribution of next point given previous point
 matrix[N,U] meen;    // fitted value for interest rate at each point
 
 //Now compute C and D functions
 
 h = (kaprn^2 + 2*beta)^0.5;
 for (t in 1:U) tau[t] = ts[t];
 for (t in 1:U) {W[t] = 2*h-(kaprn+h)*(1-exp(h*tau[t])); //tau[t] is tau
  D[1,t] = -2*(1-exp(h*tau[t]))/W[t]/tau[t];  //CIR D
  }
for (t in 1:U) {
     C[1,t] = -(omrn/beta/tau[t])*(2*log(2*h/W[t])+tau[t]*kaprn+tau[t]*h);
  }
  
// now compute priors for the short-rate process

 om = omrn;
 for (n in 1:N-1) {
Cmean[n] = rC[n]/exp(kap/opy) + om*(1-exp(-kap/opy))/kap; // mean for CIR process
Cvar[n]=beta/kap*(1-exp(-kap/opy))*(rC[n]/exp(kap/opy)+om*(1-exp(-kap/opy))/2)/kap;
   } //gamma mean and variance for next CIR step, used in prior for next rC 
 //where step d = 1/opy

// Now the parameters for the distribution of the fitted values

for (n in 1:N) {
  for (u in 1:U) {meen[n,u] = D[1,u]*rC[n]+C[1,u];  
}}
sigma_y=exp(log_sigma_y);//exp of uniform on reals makes a positive prior unbiased  
}

model { //priors and fitted distributions specified here
kaprn ~ gamma(0.9,625);
omrn ~ gamma(0.9,10);
beta ~ gamma(0.9,450);
kap ~ gamma(0.9,22.5); 
rC[1] ~ gamma(2*om/beta,2*kap/beta);  // long term distribution of short rate 
//used for first one, then incremental changes:
for (n in 2:N)  rC[n] ~ gamma(Cmean[n-1]^2/Cvar[n-1],Cmean[n-1]/Cvar[n-1]); 
//in Stan mean of gamma is alpha/beta, etc
 
for (u in 1:U) {for (n in 1:N) y[n,u] ~ normal(meen[n,u],sigma_y);}} 

generated quantities { //output log likelihood for loo
 vector[o] log_lik;
 for (j in 1:o) {
log_lik[j] = normal_lpdf(y[z[j,1], z[j,2]] | meen[z[j,1], z[j,2]], sigma_y);
}}
