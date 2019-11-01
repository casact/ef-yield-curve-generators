/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. 
 */

//Code for Single CIR ODE Case
functions{
  real[] AB_eq(real t,    // time
           real[] AB,     // A and B
		   real[] params,
           real[] x_r,    // data (real)
           int[] x_i) {   // data (integer)
  real dABdt[2];
  real A;
  real B;
  real kpa1;
  real oga1;
  real bta1;

  kpa1 = params[1];
  oga1 = params[2];
  bta1 = params[3];
  A = AB[1];
  B = AB[2];

  dABdt[1] = - oga1 * B; //dA(tau)  - delta_0
  dABdt[2] = 1 - kpa1*B - B^2*bta1/2; //dB(tau)/dtau gamma_1 + 
  
  return { dABdt[1], dABdt[2] };
  }
}

data {
  int<lower=0> N;         //# rows 
  int<lower=0> U;         //# columns 
  matrix[N,U] y;    
  real ts[U];
  real opy;          // observations per year
}

transformed data {

  int z[U*(N-1),2];  //index for rows and columns of each observation
  int o;        //# observations in y
  o = 0;
  for (n in 2:N) {
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

transformed parameters {

 real om;
 vector[U] tau;
 matrix[1,U] C;                   //C(T)
 matrix[1,U] D;                   //D(T)
 real sigma_y;     
 vector[N-1] Cmean;
 vector[N-1] Cvar;
 matrix[N,U] meen;    // fitted value at each point
 real AB0[2];
 real paras[3]; //
 real AB[U,2];   

//Compute C and D functions

 for (t in 1:U) tau[t] = ts[t];
 
 for (j in 1:2) AB0[j] = 0;
 paras[1] = kaprn;
 paras[2] = omrn;
 paras[3] = beta;
 
 AB = integrate_ode_rk45(AB_eq, AB0, 0, ts, paras, rep_array(0.0, 0), 
                         rep_array(0, 0),1e-9,1e-8,5e3);

 for (u in 1:U) {     
  C[1,u] = -AB[u,1]/tau[u];
  D[1,u] = AB[u,2]/tau[u];
 }

// Compute short rate processes  

 om = omrn;
 for (n in 1:N-1) {
  Cmean[n] = rC[n]/exp(kap/opy) + om*(1-exp(-kap/opy))/kap; 
  // mean for CIR process
  Cvar[n]=beta/kap*(1-exp(-kap/opy))*(rC[n]/exp(kap/opy)+om*(1-exp(-kap/opy))/2)/kap;
 } 
//gamma mean and variance for next CIR step, in prior for next rC;  step d = 1/opy


// Finally moments of the fit

 for (n in 1:N) {
  for (u in 1:U) {
   meen[n,u] = D[1,u]*rC[n]+C[1,u];
  }
 }
 sigma_y = exp(log_sigma_y);     //makes prior unbiased  

}

model {
kaprn ~ gamma(0.9,625);
omrn ~ gamma(0.9,10);
beta ~ gamma(0.9,450);
kap ~ gamma(0.9,22.5); 
 rC[1] ~ gamma(2*om/beta,2*kap/beta);  // long term distribution of short rate
 for (n in 2:N)  rC[n] ~ gamma(Cmean[n-1]^2/Cvar[n-1],Cmean[n-1]/Cvar[n-1]);
 //in Stan mean of gamma is alpha/beta, etc
 for (u in 1:U) {for (n in 1:N) y[n,u] ~ normal(meen[n,u],sigma_y);}
} 

generated quantities {
 vector[o] log_lik;
for (j in 1:o) {
log_lik[j] = normal_lpdf(y[z[j,1], z[j,2]] | meen[z[j,1], z[j,2]], sigma_y);
}}
