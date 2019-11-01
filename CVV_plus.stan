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

parameters {  
 real<lower=0, upper=7> kaprn1;    
 real<lower=0, upper=9> kaprn2;    
 real<lower=0, upper=9> kaprn3;   
 real<lower=0, upper=9> omrn1;   
 real omrn2;    
 real omrn3;      
 real lb;    
 real ls2;    
 real<lower=0, upper=4>  ls3;    //diff sig3 - sig2
 real<lower=-0.95, upper=-0.75> corr;     
 real<lower=-6, upper=-1> log_sigma_y;     //log stdev of residuals 
 real<lower=0, upper=7> kap1;  
 real<lower=0, upper=9> kap2;    
 real<lower=0, upper=9> kap3;   
 real om2;    
 real om3;    
 real<lower=-6, upper=-0.5> ldel;     //log of delta_0
 real gam1;
 vector<lower=0, upper=20>[N] rC;     //CIR short rate
 matrix<lower=-25, upper=25>[N,2] rV; //Vasicek short rate for each observed time point}  
}

transformed parameters {
  real beta;
 vector[3] kap;
 vector[3] om;
 vector[3] sig;
 vector[3] kaprn;  //rn for risk-neutral
 vector[3] omrn;
 vector[U] tau;
 matrix[3,U] C;                   //C(T)
 matrix[3,U] D;                   //D(T)
 vector[U] W; 
 real del;
 real sigma_y;      
 real h;
 vector[U] corradj;   
 vector[2] Cmv;  //first CIR mean and variance   
 vector[N-1] Cmean;
 vector[N-1] Cvar;
 matrix[N-1,2] Vmean; // means after first one for Vasicek processes
 vector[2] Vm;  //first Vasicek mean 
 vector[2] Vvar;  //first Vasicek variance 
 matrix[2,2] Sig;     // covariance matrix for Vasicek
 matrix[N,U] meen;    // fitted value at each point

 //Compute C and D functions
 
 del = exp(ldel);
  beta = exp(lb);
 kaprn[1] = kaprn1;
 omrn[1] = omrn1;
 sig[1] = 1;

 kaprn[2] = kaprn2;
 omrn[2] = omrn2;
 sig[2] = exp(ls2);

 kaprn[3] = kaprn3;
 omrn[3] = omrn3;
 sig[3] = sig[2]+ls3;

 for (t in 1:U) tau[t] = ts[t];
 
 h = (kaprn[1]^2 + 2*beta)^0.5;
 //W = 1/Q for CIR
 for (t in 1:U) {
  W[t] = 2*h-(kaprn1+h)*(1-exp(h*tau[t])); //tau[t] is tau
  D[1,t] = -2*(1-exp(h*tau[t]))/W[t]/tau[t]+gam1;  //CIR
  D[2,t] = (1-exp(-kaprn2*tau[t]))/kaprn2/tau[t]; //Vasicek
  D[3,t] = (1-exp(-kaprn3*tau[t]))/kaprn3/tau[t]; //Vasicek
 }
 for (t in 1:U) {
  C[1,t] = -(omrn1/beta/tau[t])*(2*log(2*h/W[t])+tau[t]*kaprn1+tau[t]*h) + del;
  C[2,t] = sig[2]^2*(-2+2*D[2,t]+kaprn[2]*tau[t]*D[2,t]^2)/4/kaprn[2]^2+omrn[2]/kaprn[2]*(1-D[2,t]);
  C[3,t] = sig[3]^2*(-2+2*D[3,t]+kaprn[3]*tau[t]*D[3,t]^2)/4/kaprn[3]^2+omrn[3]/kaprn[3]*(1-D[3,t]);
 }
 for (t in 1:U) corradj[t] = (D[2,t]+D[3,t]-1+(exp(-tau[t]*(kaprn[2]+kaprn[3]))-1)/(kaprn[2]+kaprn[3])/tau[t])*corr*sig[2]*sig[3]/kaprn[2]/kaprn[3]; 

// Compute short rate processes  

 kap[1] = kap1;
 om[1] = omrn1;
 om[2] = om2;
 kap[2] =kap2;
 om[3] = om3;
 kap[3] = kap3; 

 Cmv[1] = om[1]/kap[1];
 Cmv[2] = om[1]/kap[1]^2*beta/2;  //variance om1/kap1^2*beta/2

 for (n in 1:N-1) {
  Cmean[n] = rC[n]/exp(kap[1]/opy) + om[1]*(1-exp(-kap[1]/opy))/kap[1]; // mean for CIR process
  Cvar[n] = beta/kap[1]*(1-exp(-kap[1]/opy))*(rC[n]/exp(kap[1]/opy) + om[1]*(1-exp(-kap[1]/opy))/2)/kap[1];
 } //gamma mean and variance for next CIR step, used in prior for next rC where step d = 1/opy

 Vm[1] = om[2]/kap[2];
 Vm[2] = om[3]/kap[3];
 Vvar[1] = sig[2]^2/2/kap[2];
 Vvar[2] = sig[3]^2/2/kap[3];

 for (n in 1:N-1){ Vmean[n,1] = rV[n,1]+(om[2]-kap[2]*rV[n,1])/opy;
                   Vmean[n,2] = rV[n,2]+(om[3]-kap[3]*rV[n,2])/opy;} //mean for rV[n+1,i]
 for (i in 1:2) Sig[i,i] = sig[i+1]^2/opy;
 Sig[1,2] = corr*sig[3]*sig[2]/opy;
 Sig[2,1] = Sig[1,2];

// Finally moments of the fit
 for (n in 1:N) {
  for (u in 1:U) {
   meen[n,u] = D[1,u]*rC[n]+C[1,u]+D[2,u]*rV[n,1]+C[2,u]+D[3,u]*rV[n,2]+C[3,u]+corradj[u];
  }
 }
 sigma_y = exp(log_sigma_y);     //makes prior unbiased  
 
//Priors not specified are uniform over their defined ranges, like non-negative reals for sigma_y; transformed parameters are simulated but don't have priors.
}

model {
kaprn1 ~ gamma(160,250);
kaprn2 ~ gamma(1400,38000);
kaprn3 ~ gamma(400,6600);
omrn1 ~ gamma(20,6);
omrn2 ~ normal(9.5,0.2);
omrn3 ~ normal(-9,0.35);
lb ~ normal(-0.5,0.1);
ls2 ~ normal(0.12,0.05);
ls3 ~ normal(0.13,0.01);
kap1 ~ gamma(195,140);
kap2 ~ gamma(1000,1000);
kap3 ~ gamma(20,45);
om2 ~ normal(3,1);
om3 ~ normal(-1.6,0.5);
gam1 ~ double_exponential(0,0.01);
 rC[1] ~ gamma(2*om[1]/beta,2*kap[1]/beta);  // long term distribution of short rate
 for (n in 2:N)  rC[n] ~ gamma(Cmean[n-1]^2/Cvar[n-1],Cmean[n-1]/Cvar[n-1]); //in Stan mean of gamma is alpha/beta, etc
 rV[1,1] ~ normal(Vm[1],Vvar[1]^0.5); 
 rV[1,2] ~ normal(Vm[2],Vvar[2]^0.5); 
for (n in 2:N) {rV[n,] ~ multi_normal(Vmean[n-1,],Sig);}
 for (u in 1:U) {for (n in 1:N) y[n,u] ~ normal(meen[n,u],sigma_y);}
} 

generated quantities {
 vector[o] log_lik;
for (j in 1:o) log_lik[j] = normal_lpdf(y[z[j,1], z[j,2]] | meen[z[j,1], z[j,2]], sigma_y);
}
