

// ###################  SIR  ODE Function  
  
  functions {
    
    real[] SIR(real t , real[] y , real[] parms , real[] rdata , int[] idata){
      
      real beta = parms[1];
      real gamma = parms[2];
      real rho = parms[3];
      real dydt[3];
      
      dydt[1] = -beta*y[1]*y[2];
      dydt[2] = beta*y[1]*y[2]-gamma*y[2];
      dydt[3] = gamma*y[2];
      
      
      return dydt;
    }
  }


// ################### Data Module 
  
  data {
    
    int<lower=0> k;   //number of infected individuals 
    real<lower=0> t0;     // initial time 
    real<lower=0> ti[k];  // times of infection 
    real<lower=0> Ti[k];    // times of recovery 
    int<lower=0,upper=1> event[k];   // censoring variable;  default =1 (no censoring) 
    
  }

transformed data{
  
  real x_r[0];
  int x_i[0];
  
}

parameters {
  real<lower=0.15> beta;               
  real<lower=0, upper=beta> gamma;               
  real<lower=0, upper=1> rho;   
}

transformed parameters{
  real R0 = beta/gamma;
  real N=1/rho;
  real c = beta*rho;
}


// ##################### Model Defining Module  
  
  model {
    
    real temp[k,3];
    real parms[3];
    real init[3];
    
    
    
    parms[1] = beta; 	//rate of infection 
    parms[2]=gamma;  	// rate of recovery
    parms[3]= rho;	// initial infection prop. 
    
    init[1] = 1;
    init[2] = rho; 
    init[3] = 0;
    
    // ## likelihood function based on SIR ODE 
      
      temp = integrate_ode_rk45(SIR,init,t0,ti,parms,x_r,x_i,1.0E-7,1.0E-7,1.0E7); //-6
    
    for(i in 1:k){
      target += log(temp[i,1])+log(temp[i,2])+event[i]*log(gamma)-gamma*(Ti[i]-ti[i]);
      
    }
    
    target += k*(log(beta)-log(1-temp[k,1])) ;
    target += gamma_lpdf(beta|.002,.002)+gamma_lpdf(gamma|.002,.002);  
    
    // ## non-informative priors 
      
  }
  
  
