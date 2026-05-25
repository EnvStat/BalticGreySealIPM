
functions { // function for computing the expected composition of the hunting bag 
  array[] real dH_dtau(real tau,array[] real H, array[] real vars, data array[] real d_r, data array[]  int d_i) {
    array[1] real dH;
    real n0=vars[1];
    real k=vars[2];
    real E_1=vars[3];
    real E_2=vars[4];
    real mu=vars[5];
    dH[1] = n0*exp(-(E_1+E_2)*(k*tau-tau^2/2) - mu*tau)*E_1*(k-tau);
    return dH;
  }
}

data {
  //years and age classes
  int<lower=0> t; // years to run the modle
  int<lower=0> a;

  //initial population structure
  vector[2*a] N0;
  
  //Prior covariance matrix for hunting selectivities (called biases in the code)
  matrix[2*a,2*a] L_bias;

  //Herring
  vector[t + 1]  h_BP_GoF;
  vector[t + 1]  h_GoB;

  //Hunting quota Sweden
  array[t] int <lower=0> Q_sw;
  array[t] int <lower=0> Q_fi;
  
  //Hunting bags Sw&Fin
  array[t] real <lower=0> y_hb_sw;
  array[t] real <lower=0> y_hb_fi;
  
  //Hunting samples Sw&Fin
  array[2*a,t] int <lower=0> y_hs_sw;
  array[2*a,t] int <lower=0> y_hs_fi;

  //Bycatch samples 
  array[2*a,t] int <lower=0> y_bc;
  
  //Aerial observations
  array[t-1] int<lower=0> y;
  
  //Pregnancy observatiuons
  array[t] int<lower=0> y_pr;
  array[t] int<lower=0> y_pr_tot;
  
  //Reproductive signs
  array [4, t] int<lower=0> z_fi;
}

transformed data{
  
  //Aging matrix
  matrix[2*a,2*a] A;
  for (i in 1:2*a){
    for (j in 1:2*a){
      if(i-j==1){
        A[i,j]=1;
      }
      else{
        A[i,j]=0;
      }
    }
  }
  A[a,a]=1;
  A[a+1,a]=0;
  A[2*a,2*a]=1;
  
  
  real tau_s=1.0/36; // Time between mating and pregnancy sample collection
  real tau_h_0=1.5/12;  //Time between seals giving birth and start of the hunting season on 20th of April
  real tau_h=8.0/12;  // Time between seals giving birth and end of hunting season, this is k in the model
  real tau_h_m=5.5/12; //Length of the hunting season
  
  array[0] real empty_real_array;
  array[0] int empty_int_array;

}

parameters {
  
  // Biases (hunting selectivity)
  vector[2*a] g_sw_sc;
  vector[2*a] g_fi_sc;
  vector[2*a] g_bc_sc;
  
  
  real<lower=0> sigma_h_sw;
  real<lower=0> sigma_h_fi;
  vector<lower=0> [t] epsilon_h_sw;
  vector<lower=0> [t] epsilon_h_fi;


  //Carrying capacity
  real<lower=0> K;

  //Natural mortality
  real<lower=0, upper=1> phi_a_sc; 
  real<lower=0, upper=1> phi_sc;
  real v0;
  real v5;
  real <lower=0, upper=1> c;

  //Birth rate
  real<lower=0, upper=1> b0max;
  real<lower=0, upper=1> b0min_sc;

  real alpha_sc;

  real<lower=0, upper=1> theta0_sc;
  real beta;
  real<lower=0, upper=1> w;


  //State process
  real n0;
  vector[t+1] epsilon_birth;
  vector[t+1] epsilon_sex;
  matrix[3*2*a,t] u;

  //Aerial observation
  real<lower=0, upper=1> mu;
  real<lower=0> r; 

  //Reproductive signs observation
  real<lower=0, upper=1> kappa;
  real<lower=0, upper=1> pi_c_mean;
  real<lower=0, upper=1> pi_s_mean;
  vector<lower=0>[t] epsilon_rep_c;
  vector<lower=0>[t] epsilon_rep_s;
  real<lower=0> sigma_rep_c;
  real<lower=0> sigma_rep_s;
}

transformed parameters{

  //Biases (hunting selectivity)
  vector[2*a] g_sw=L_bias*g_sw_sc;
  vector[2*a] g_fi=L_bias*g_fi_sc;
  vector[2*a] g_bc=L_bias*g_bc_sc;

  array[5] real vars;
  vars[2]=tau_h;
  
  array[1] real init_state;
  init_state[1]=0.0;
  array[1] real times;
  times[1]=tau_h;

  
  
  vector<lower=0>[2*a] E_sw;
  vector<lower=0>[2*a] E_fi;
  
  vector<lower=0> [2*a] N_h_sw;
  vector<lower=0> [2*a] N_h_fi;


  
  
  //Mortality
  vector  <lower=0>[2*a] mu_m;
  real <lower=0, upper=1> phi_a=0.8+0.2*phi_a_sc;
  real <lower=0, upper=1> phi_pup;


  //Hunting 

  //Sweden
  
  matrix[2*a,t] H_sw;
  vector<lower=0>[t] H_sw_tot;
  
  
  matrix[2*a,t] D;
  vector<lower=0>[t] D_tot;
  
  //Finland
  matrix[2*a,t] H_fi;
  vector<lower=0>[t]  H_fi_tot; 
  
  
  //Birth rate
  real alpha=alpha_sc*beta;

  vector<lower=0,upper=1>[t+1]  b0;

  array[t+1] real<lower=0,upper=1>  b;
  real <lower=0, upper=0.8> bK;
  real <lower=0> theta0;
  real <lower=0, upper=1>  b0min=b0max*b0min_sc;
  real <lower=0, upper=1> b0_av=b0min+ (b0max-b0min)*inv_logit(alpha);


  real theta1;
  
  //Preganancy rate
  array[t] real<lower=0, upper=1> p;
  
  //State process
  matrix[2*a,t+1] N;
  vector[2*a] S_diag;
  array[t+1] real<lower=0> N_temp;
  array [t+1] real N_tot; 
  array[t] matrix [2*a, 4] rho;
    
  //Transformation matrix
  matrix [2*a,2*a] M_S;
  matrix [2*a,2*a] M_D;
  matrix [2*a,2*a] M_H_sw;
  matrix [2*a,2*a] M_H_fi;

  matrix [4*2*a,2*a] M;
  
  
  //Helper matrices
  matrix[8*a, t] U;
  matrix[2*a, 4] U_stacked;
  matrix [2*a,3] u_stacked;
    
  matrix [2*a, 4] eta_p;
  row_vector[4] eta_p_adj;
    
  
  row_vector[3] Mean;
  matrix [3,3] Sigma;
  matrix[3,3] L;
    
  row_vector[4] allocation;

  //Reproductive sign probabilities
  matrix[4,t] gamma;
  
 

  vector<lower=0,upper=1>[t] pi_c=pi_c_mean*exp(-epsilon_rep_c*sigma_rep_c);
  vector<lower=0,upper=1>[t] pi_s=pi_s_mean*exp(-epsilon_rep_s*sigma_rep_s);
  
  

  // Birth rate

  b0 = b0min+ (b0max-b0min)*inv_logit(alpha+beta*(w*h_BP_GoF+(1-w)*h_GoB));
    
  theta0 = -log(b0max+(1-b0max)*theta0_sc);
 
  
  //Natural mortality
  phi_pup=phi_a*phi_sc;
  mu_m[1]=-log(phi_pup);
  mu_m[6]=-log(phi_a);
  
  for (j in 2:5){
      mu_m[j]=exp(log(mu_m[1])+((j-1)*1.0/5)^c*(log(mu_m[6])-log(mu_m[1])));
  }
  mu_m[7]=exp(log(mu_m[1])+v0);
  mu_m[12]=exp(log(mu_m[6])+v5);

  for (j in 2:5){
      mu_m[6+j]=exp(log(mu_m[7])+((j-1)*1.0/5)^c*(log(mu_m[12])-log(mu_m[7])));
  }
  S_diag=exp(-mu_m);

  //Birth rate at carrying capacity
  bK=2*(1-phi_a)/exp(sum(-mu_m[1:5]));
  theta1=log(1-log(bK/b0_av)/(theta0))/K;

  for (i in 1:t){
    if (i==1){ //Year 2002 april
    
      // Birth rate
      b[i]=b0[i]*exp(-theta0*(exp(theta1*sum(N0))-1));
      
      //State processes
      N[,i]=N0;
      
      for (k in 1:20){
        N[,i]=A*diag_matrix(S_diag)*N[,i];
        N[1,i]=b[i]/2*N[6,i];
        N[7,i]=b[i]/2*N[6,i];
      }
      
      N[,i]=N[,i]*n0/sum(N[,i]);
      N_tot[i]=sum(N[,i]);
      N_temp[i]=N[1,i]; 
    }
    else{
      //Birth rate
      b[i]=b0[i]*exp(-theta0*(exp(theta1*N_tot[i-1])-1));

      //Pregnancy rate
      p[i-1]= b0[i]*exp(theta0*(1-tau_s*exp(theta1*N_tot[i-1])));


      //State process
      
      N[,i]=A*U_stacked[,1];
      N_temp[i]=N[a,i]*b[i]+sqrt(N[a,i]*b[i]*(1-b[i]))*epsilon_birth[i];
      N[1,i]=N_temp[i]/2+sqrt(N_temp[i]/4)*epsilon_sex[i];
      N[a+1,i]=N_temp[i]-N[1,i];

      N_tot[i]=sum(N[,i]);
      
    }
    
    //Expected number of hunted seals in each demographic group
   
    E_sw=(Q_sw[i]*exp(-epsilon_h_sw[i]*sigma_h_sw)*2/(tau_h^2))*exp(g_sw)/(exp(g_sw)'*N[,i]);
    E_fi=(Q_fi[i]*exp(-epsilon_h_fi[i]*sigma_h_fi)*2/(tau_h^2))*exp(g_fi)/(exp(g_fi)'*N[,i]);
  
    for (x in 1:(2*a)){
       vars[1]=N[x,i];
       vars[5]=mu_m[x];

      vars[3]=E_sw[x];
      vars[4]=E_fi[x];
  
      N_h_sw[x] = integrate_ode_rk45(dH_dtau, init_state,0.0, times,
    vars,empty_real_array, empty_int_array)[1,1];
    }
    
        
    for (x in 1:(2*a)){
       vars[1]=N[x,i];
       vars[5]=mu_m[x];

      vars[4]=E_sw[x];
      vars[3]=E_fi[x];
      N_h_fi[x] = integrate_ode_rk45(dH_dtau, init_state,0.0, times,
    vars,empty_real_array, empty_int_array)[1,1];
    }
    
    
    
    rho[i,,3]=N_h_sw ./ N[,i];
    rho[i,,4]=N_h_fi./ N[,i];
    rho[i,,1]=exp(-(E_sw+E_fi)*tau_h^2/2) .* S_diag;
    rho[i,,2]=1-rho[i,,1]-rho[i,,3]-rho[i,,4];
   

    // Transition matrix
    M_S=diag_matrix(rho[i,,1]);
    M_D=diag_matrix(rho[i,,2]);
    M_H_sw=diag_matrix(rho[i,,3]);
    M_H_fi=diag_matrix(rho[i,,4]);

    M=append_row(M_S, append_row(M_D, append_row(M_H_sw, M_H_fi)));
    
    //Multinomial transition modelled as logit-N
    
    eta_p=to_matrix(M*N[,i], 2*a, 4);
    

    u_stacked=to_matrix(u[,i],2*a, 3);
    for (j in 1:(2*a)){
      
      eta_p_adj=eta_p[j,]*(1+1/min(eta_p[j,1:4]));
      Mean=digamma(eta_p_adj[2:4])-digamma(eta_p_adj[1]);
      Sigma=rep_matrix(trigamma(eta_p_adj[1]),3,3)+diag_matrix(trigamma(eta_p_adj[2:4])');
 
 
      L=cholesky_decompose(Sigma);
      allocation=softmax(append_col(0, Mean+u_stacked[j,]*L')')';
      U_stacked[j,]=allocation*N[j,i];
    }
    U[,i]=to_vector(U_stacked);
    
    // Actual number of hunted and bycaught seals
    H_sw[,i]=U_stacked[,3];
    H_sw_tot[i]=sum(H_sw[,i]);
    
    H_fi[,i]=U_stacked[,4];
    H_fi_tot[i]=sum(H_fi[,i]);
    
    D[,i]=U_stacked[,2];
    D_tot[i]=sum(D[,i]);

    //Observation probabilities of reproductive signs
    gamma[2,i]=b[i]*pi_s[i]*(1-pi_c[i]); //gamma_10 PS yes, CA no
    gamma[3,i]=b[i]*(1-pi_s[i])*pi_c[i]+(1-b[i])*kappa*pi_c[i];
    gamma[4,i]=b[i]*pi_s[i]*pi_c[i];
    gamma[1,i]=1-sum(gamma[2:4,i]);
  }
    
    // Calculate birth and pregnancy rates and population size for 2024
  
    //Birth rate
    b[t+1]=b0[t+1]*exp(-theta0*(exp(theta1*N_tot[t])-1));

    //Pregnancy rate
    p[t]=b0[t+1]*exp(theta0*(1-tau_s*exp(theta1*N_tot[t])));


    N[,t+1]=A*U_stacked[,1];
    N_temp[t+1]=N[a,t+1]*b[t+1]+sqrt(N[a,t+1]*b[t+1]*(1-b[t+1]))*epsilon_birth[t+1];
    N[1,t+1]=N_temp[t+1]/2+sqrt(N_temp[t+1]/4)*epsilon_sex[t+1];
    N[a+1,t+1]=N_temp[t+1]-N[1,t+1];
    N_tot[t+1]=sum(N[,t+1]);
}

model {
  //Initial population size
  n0 ~ lognormal (9.8, 0.1);

 //Natural mortality
  phi_a_sc ~ uniform(0,1);
  phi_sc ~ uniform (0,1);
  c ~uniform (0,1);
  v0 ~ cauchy(0,0.2);
  v5 ~ cauchy(0.88,0.2);
  
  //Hunting and bycatch bias (selectivity)
  
  g_sw_sc ~ normal(0,0.5);
  g_fi_sc ~ normal(0,0.5);
  g_bc_sc ~ normal(0,0.5);
  
  //Hunting effort sd
  sigma_h_sw ~cauchy(0,0.1);
  sigma_h_fi ~ cauchy(0,0.1);
  
  //Birth rate
  b0max ~ uniform (0,1);
  b0min_sc ~ uniform (0,1);
  
  alpha_sc ~ normal(0,4);
  beta ~ normal(0,3);
  w ~ uniform(0,1);

  theta0_sc ~ uniform(0,1);
  
  //Carrying capacity
  K ~ lognormal(11.3,0.3);


  //Observation of aerial survey  
  mu ~ beta(32,9);
  r ~ lognormal (5.3, 1);

  //Observtion of reproductive signs
  kappa ~ uniform(0,1);
  pi_s_mean ~ uniform(0,1);
  pi_c_mean ~ uniform(0,1);
  
  sigma_rep_c ~ normal (0,0.1);
  sigma_rep_s ~ normal (0,0.1);
  
  //standard normals for stochasticity
  
  epsilon_h_sw ~ normal(0,1);
  epsilon_h_fi ~normal(0,1);

  epsilon_rep_c ~normal(0,1);
  epsilon_rep_s ~normal(0,1);

  //Stochasicity for birth process
  epsilon_birth ~ normal(0,1);
  epsilon_sex ~ normal(0,1);
  
  //Stochasticity for state transitions
  for (i in 1:t){
    u[,i] ~ normal(0,1);
  }
  
  
  ////////////////////OBSERVATION MODEL ////////////////////////

  for ( i in 1:t){
  //Aerial observations
    if(i<=22){
      if (y[i]>10){
        y[i] ~ neg_binomial_2(mu*N_tot[i+1],r);
      }
    }
  
    //Hunting bag
    y_hb_sw[i] ~ normal(H_sw_tot[i], 0.05*H_sw_tot[i]);
    y_hb_fi[i] ~ normal(H_fi_tot[i], 0.05*H_fi_tot[i]);

    
  // Hunting samples
   y_hs_sw[,i] ~ multinomial(H_sw[,i]/H_sw_tot[i]);

   if (i<22){
     y_hs_fi[,i] ~ multinomial(H_fi[,i]/H_fi_tot[i]);
   }

   //By-catch samples
   y_bc[,i] ~ multinomial(exp(g_bc).*D[,i]/(exp(g_bc)'*D[,i]));
    
    
    //Pregnancy
    y_pr[i] ~ binomial( y_pr_tot[i], p[i]);
    
    //Reproductive signs
    if (i<=22){
      z_fi[,i] ~ multinomial(gamma[,i]);
    }
  }
}

generated quantities{
  array[t] real<lower=0> y_hb_sw_sim;
  array[t] real<lower=0> y_hb_fi_sim;

  //Aerial observations
  array[t+1] int<lower=0> y_sim;

  //Pregnancy observatiuons
  array[t] int<lower=0> y_pr_sim;


  //Reproductive signs
  array[4, t ] int<lower=0> z_sim;
  

  y_sim=neg_binomial_2_rng(mu*to_vector(N_tot), r);

   for ( i in 1:t){
     y_pr_sim[i]=binomial_rng( y_pr_tot[i], p[i]);
     if (sum(z_fi[,i])>0){
      z_sim[,i]= multinomial_rng(gamma[,i], sum(z_fi[,i]));
     }
     else{
       z_sim[,i]=rep_array(0,4);
     }

   }
   y_hb_sw_sim = normal_rng(H_sw_tot, 0.05*to_vector(H_sw_tot));
   y_hb_fi_sim = normal_rng(H_fi_tot, 0.05*to_vector(H_fi_tot));

}


