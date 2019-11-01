/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. 
 */

 functions{
  real[] AB_eq(real t,        // time
           real[] AB,      // A and B1 to B3
		   real[] params,
           real[] x_r,    // data (real)
           int[] x_i) {   // data (integer)
  real dABdt[6];
  real A[3];
  real B[3];
  real kpa1;
  real kpa2;
  real kpa3;
  real oga1;
  real oga2;
  real oga3;
  real bta1;
  real bta2;
  real bta3;
  real sigma2;
  real sigma3;
  real sigma32;
  real kpa21;
  real kpa23;
  real kpa31;
  real kpa32;

  kpa1 = params[1];
  oga1 = params[2];
  bta1 = params[3];
  kpa2 = params[4];
  oga2 = params[5];
  bta2 = params[6];
  sigma2 = params[7];
  kpa3 = params[8];
  oga3 = params[9];
  bta3 = params[10];
  sigma3 = params[11];
  sigma32 = params[12];
  kpa21 = params[13];
  kpa23 = params[14];
  kpa31 = params[15];
  kpa32 = params[16];
  A[2] = AB[2];
  A[3] = AB[3];
  B[1] = AB[4];
  B[2] = AB[5];
  B[3] = AB[6];

  dABdt[1] = - oga1 * B[1]; //dA(tau)  - delta_0
  dABdt[2] = - oga2 * B[2] + 0.5 * (sigma2 * B[2] + sigma32 * B[3])^2; //dA(tau)
  dABdt[3] = - oga3 * B[3] + 0.5 * (sigma3 * B[3])^2 ; //dA(tau)
  dABdt[4] = 1 - kpa1*B[1] - kpa21*B[2] - kpa31*B[3] - B[1]^2*bta1/2; //dB(tau)/dtau gamma_1 + 
  dABdt[5] = 1 - kpa2*B[2] - kpa32*B[3] - (sigma2*B[2]+sigma32*B[3])^2*bta2/2; //dB(tau)/dtau
  dABdt[6] = 1 - kpa3*B[3] - kpa23*B[2] - sigma3^2*B[3]^2*bta3/2; //dB(tau)/dtau
  
  return { dABdt[1], dABdt[2], dABdt[3], dABdt[4] , dABdt[5], dABdt[6] };
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

  int z[U*N,2];  //index for rows and columns of each observation
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

parameters	{												
 real<lower=-0, upper=1> kaprn1;  //CIR 
 real<lower=-0, upper=10> kaprn2;  //Vasicek
 real<lower=-0, upper=10> kaprn3;  //Vasicek
 real kaprn21;    
 real kaprn23;     
 real kaprn31;    
 real kaprn32;    
 real<lower=-0, upper=10> omrn1;      //CIR
 real omrn2;      //Vasicek 
 real omrn3;      //Vasicek
 real lb1;   
 real lb2;     
 real lb3;   
 real ls2;    //sigma for 1st Vasicek
 real ls3;    //sigma for 2nd Vasicek
 real<lower=-0.99, upper=0> corr;       
 real<lower=-6, upper=-1> log_sigma_y;//log stdev of residuals 
 real<lower=-0, upper=1> kaprw1;
 real<lower=-0, upper=10> kaprw2;
 real<lower=-0, upper=10> kaprw3;
 real kap21;
 real kap23;
 real kap31;
 real kap32;
 real omrw2;
 real omrw3;
 real<lower=-3, upper=-0.4> ldel; //log of delta_0
 real gam1;
 vector<lower=0, upper=20>[N] rC; //CIR short rate
 matrix<lower=-25, upper=25>[N,2] rV;  //Vasicek short rates 
}													

transformed parameters {
 vector[3] beta;
 vector[3] kap;
 vector[3] om;
 vector[3] sig;
 vector[3] kaprn;
 vector[3] omrn;
 vector[U] tau;
 matrix[3,U] C;                 
 matrix[3,U] D;                 
 real sig_32;
 real del;
 real sigma_y;      
 vector[2] Cmv;  //first CIR mean and variance   
 vector[N-1] Cmean;
 vector[N-1] Cvar;
 matrix[N-1,2] Vmean; // means after first one for Vasiceks
vector[2] Vm;     //first Vasicek mean 
 vector[2] Vvar;  //first Vasicek variance 
matrix[N,U] meen; // fitted value at each point
 real AB0[6];
 real paras[16]; 
 real AB[U,6];   
 real w2;
 real w3;
 real kd;

//Compute C and D functions
 
 del = exp(ldel);
 beta[1] = exp(lb1);
 kaprn[1] = kaprn1;
 omrn[1] = omrn1;
 sig[1] = 1;

 beta[2] = exp(lb2);
 kaprn[2] = kaprn2;
 omrn[2] = omrn2;
 sig[2] = exp(ls2);
 sig_32 = corr*sig[2];

 beta[3] = exp(lb3);
 kaprn[3] = kaprn3;
 omrn[3] = omrn3;
 sig[3] = (1-corr^2)^0.5*exp(ls3);
 
 for (t in 1:U) tau[t] = ts[t];

 for (j in 1:6) AB0[j] = 0;
 paras[1] = kaprn1;
 paras[2] = omrn1;
 paras[3] = beta[1];
 paras[4] = kaprn2;
 paras[5] = omrn2;
 paras[6] = beta[2];
 paras[7] = sig[2];
 paras[8] = kaprn3;
 paras[9] = omrn3;
 paras[10] = beta[3];
 paras[11] = sig[3];
 paras[12] = sig_32;
 paras[13] = kaprn21;
 paras[14] = kaprn23;
 paras[15] = kaprn31;
 paras[16] = kaprn32;

 AB = integrate_ode_rk45(AB_eq, AB0, 0, ts, paras, rep_array(0.0, 0), rep_array(0, 0),1e-7,1e-6,5e3);

 for (u in 1:U) {     
  C[1,u] = -AB[u,1]/tau[u] + del;
  C[2,u] = -AB[u,2]/tau[u];  
  C[3,u] = -AB[u,3]/tau[u];  
  D[1,u] = AB[u,4]/tau[u] + gam1;
  D[2,u] = AB[u,5]/tau[u];
  D[3,u] = AB[u,6]/tau[u];
 }

// Compute short rate processes  

 kap[1] = kaprw1; 
 om[1] = omrn1;

 om[2] = omrw2; 
 kap[2] = kaprw2; 

 om[3] = omrw3; 
 kap[3] = kaprw3;

 w2 = om[2] - kap21*om[1]/kap[1];
 w3 = om[3] - kap31*om[1]/kap[1];
 kd = kap[2]*kap[3] - kap23*kap32;

 Cmv[1] = om[1]/kap[1];
 Cmv[2] = om[1]/kap[1]^2*beta[1]/2;  
 for (n in 1:N-1) {
  Cmean[n] = rC[n]/exp(kap[1]/opy) + om[1]*(1-exp(-kap[1]/opy))/kap[1]; // mean for CIR process
  Cvar[n] = beta[1]/kap[1]*(1-exp(-kap[1]/opy))*(rC[n]/exp(kap[1]/opy) + om[1]*(1-exp(-kap[1]/opy))/2)/kap[1];
 } //gamma mean and variance for next CIR step, used in prior for next rC where step d = 1/opy

 Vm[1] = (kap[3]*w2 - kap23*w3)/kd;
 Vm[2] = (kap[2]*w3 - kap32*w2)/kd;
 Vvar[1] = sig[2]^2*(beta[2]*om[1]/kap[1]+1)/2/kap[2];
 Vvar[2] = (sig[3]^2*(beta[3]*om[1]/kap[1]+1)+sig_32^2*(beta[2]*om[1]/kap[1]+1))/2/kap[3];
 for (n in 1:N-1){
Vmean[n,1] = rV[n,1]+(om[2]-kap21*rC[n] -kap[2]*rV[n,1] - kap23*rV[n,2])/opy;
  Vmean[n,2] = rV[n,2]+(om[3]-kap31*rC[n] -kap32*rV[n,1] - kap[3]*rV[n,2])/opy;
 } 

// Finally moments of the fit
 for (n in 1:N) {
  for (u in 1:U) {
   meen[n,u] = D[1,u]*rC[n]+C[1,u]+D[2,u]*rV[n,1]+C[2,u]+D[3,u]*rV[n,2]+C[3,u];
  }
 }
 sigma_y = exp(log_sigma_y);     

}

model { //alternate priors gave better fit but has narrow ranges and fit gets very bad outside of them
//	kaprn1 ~ gamma(0.9,90);
//	kaprn2 ~ normal(0.52,0.007);
//	kaprn3 ~ normal(0.32,0.008);
//	kaprn21 ~ normal(0.624,0.0025);
//	kaprn23 ~ normal(1.485,0.006);
//	kaprn31 ~ normal(0.12,0.0025);
//	kaprn32 ~ normal(-0.0465,0.005);
//	omrn1 ~ normal(0.16,0.009);
//	omrn2 ~ normal(1.9,0.0025);
//	omrn3 ~ normal(2,0.025);
//	lb1 ~ normal(-7.8, 0.08);
//	lb2 ~ normal(0.68,0.025);
//	lb3 ~ normal(2,0.0125);
//	ls2 ~ normal(-1.85,0.0125);
//	ls3 ~ normal(-0.54,0.03);
//	kaprw1 ~ gamma(0.9, 45);
//	kaprw2 ~ normal(0.87, 0.0175);
//	kaprw3 ~ normal(0.375, 0.0075);
//	kap21 ~ normal(-0.043, 0.00125);
//	kap23 ~ normal(0.48, 0.00125);
//	kap31 ~ normal(0.875, 0.0125);
//	kap32 ~ normal(0.65, 0.01);
//	omrw2 ~ normal(-6.3, 0.0125);
//	omrw3 ~ normal(2.75, 0.0125);
kaprn1 ~ gamma(12,3600);
kaprn2 ~ gamma(250,1500);
kaprn3 ~ gamma(675,5000);
kaprn21 ~ normal(0.12,0.005);
kaprn23 ~ normal(0.2,0.005);
kaprn31 ~ normal(-0.63,0.02);
kaprn32 ~ normal(-0.33,0.01);
omrn1 ~ gamma(1000,9000);
omrn2 ~ normal(-2.25,0.05);
omrn3 ~ normal(2.7,0.05);
lb1 ~ normal(-7.5,0.05);
lb2 ~ normal(0.9,0.05);
lb3 ~ normal(0.8,0.1);
ls2 ~ normal(-1.5,0.05);
ls3 ~ normal(-1.75,0.1);
kaprw1 ~ gamma(280, 12000);
kaprw2 ~ gamma(40,70);
kaprw3 ~ gamma(250,300);
kap21 ~ normal(0.2, 0.05);
kap23 ~ normal(-1, 0.7);
kap31 ~ normal(0, 0.04);
kap32 ~ normal(-0.3, 0.04);
omrw2 ~ normal(-7, 0.75);
omrw3 ~ normal(3.8, 0.25);
 gam1 ~ double_exponential(0,0.01);
 rC[1] ~ gamma(2*om[1]/beta[1],2*kap[1]/beta[1]);  // long term distribution of short rate
 for (n in 2:N)  rC[n] ~ gamma(Cmean[n-1]^2/Cvar[n-1],Cmean[n-1]/Cvar[n-1]); //in Stan mean of gamma is alpha/beta, etc
 rV[1,1] ~ normal(Vm[1],Vvar[1]^0.5); 
 rV[1,2] ~ normal(Vm[2],Vvar[2]^0.5); 
 for (n in 2:N) {
  matrix[2,2] Sig;
  Sig[1,1] = sig[2]^2*(1 + beta[2] * rC[n-1])/opy;
  Sig[2,2] = sig_32^2*(1 + beta[2] * rC[n-1])/opy + sig[3]^2*(1 + beta[3] * rC[n-1])/opy; 
 Sig[1,2] = sig_32*sig[2]*(1+beta[2]*rC[n-1])/opy;
  Sig[2,1] = Sig[1,2];
  rV[n,] ~ multi_normal(Vmean[n-1,],Sig);
 }
 for (u in 1:U) {for (n in 1:N) y[n,u] ~ normal(meen[n,u],sigma_y);}
} 

generated quantities {
 vector[o] log_lik;
for (j in 1:o) log_lik[j] = normal_lpdf(y[z[j,1], z[j,2]] | meen[z[j,1], z[j,2]], sigma_y);
}
