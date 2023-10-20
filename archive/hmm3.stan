

functions {

  matrix generate_P_matrix(
    real t, 
    vector q0, vector q1, vector q2,
    vector beta_01, vector beta_02, 
    vector beta_10, vector beta_12, 
    vector beta_20, vector beta_21,
    vector covs,
    int debug) {

    matrix[3,3] Q;
    matrix[3,3] P;

    Q[1,2] = q0[1] * exp(dot_product(beta_01, covs));
    Q[1,3] = q0[2] * exp(dot_product(beta_02, covs));
    Q[1,1] = 0 - Q[1,2] - Q[1,3];

    Q[2,1] = q1[1] * exp(dot_product(beta_10, covs));
    Q[2,3] = q1[2] * exp(dot_product(beta_12, covs));
    Q[2,2] = 0 - Q[2,1] - Q[2,3];
    
    Q[3,1] = q2[1] * exp(dot_product(beta_20, covs));
    Q[3,2] = q2[2] * exp(dot_product(beta_21, covs));
    Q[3,3] = 0 - Q[3,1] - Q[3,2];

    if (debug == 1) {
    print("Q matrix: ", Q);
    }

    P = matrix_exp(Q * t);
    return(P);
  }

  matrix generate_E_matrix(int state) {
    matrix[3,3] E;
  
    E = rep_matrix(0, 3, 3);

    if (state == 0) {
      E[1,1] = 1;
    } else if (state == 1) {
      E[2,2] = 1;
    } else if (state == 2) {
      E[3,3] = 1;
    } else if (state == 3) {
      E[2,2] = 1;
      E[3,3] = 1;
    } else if (state == 4) { 
      E[1,1] = 1;
      E[2,2] = 1;
      E[3,3] = 1;
    }

    return(E);
  }

}

data {
  // ... declarations ...
  int debug;
  int N; // number of participants
  int n; // total number of observations
  array[N] int n_obs; // number of observation-pairs per participant
  array[n] real t_start; 
  array[n] real t_end;
  // states 0 = no org
  // 1 = S
  // 2 = R
  // 3 = 1 or 2
  // 4 = 1 or 2 or 3 (ie censored)
  array[n] int state_start;
  array[n] int state_end;
  // covariate values for segment
  array[n] vector[1] covs;

}

transformed data {
   // ... declarations ... statements ...
}

parameters {
   // state transition
   //vector<lower=0>[2] q0; // state transition instantaneous probs (q01, q02)
   //vector<lower=0>[2] q1; // state transition instantaneous probs (q10, q12)
   //vector<lower=0>[2] q2; // state transition instantaneous probs (q20, q21)
   
   real<lower=0> q01; // fix same as q02
   real<lower=0> q10; // fix same as q20
   real<lower=0> q12;
   real<lower=0> q21;
   vector[1] beta_01; // covariates on state transitions beta_ab from a->b
   vector[1] beta_10; //
   //vector[1] beta_02;
   //vector[1] beta_20;
   vector[1] beta_12;
   vector[1] beta_21;
   array[N] simplex[3] p0; // start state probabilities
   //vector<lower=0>[2] beta_parameters;
}

transformed parameters {
   vector<lower=0>[2] q0; // state transition instantaneous probs (q01, q02)
   vector<lower=0>[2] q1; // state transition instantaneous probs (q10, q12)
   vector<lower=0>[2] q2; // state transition instantaneous probs (q20, q21)
   vector[1] beta_02;
   vector[1] beta_20;

   q0 = [q01, q01]';
   q1 = [q10, q12]';
   q2 = [q10, q21]';

   beta_02 = beta_01;
   beta_20 = beta_10;
}

model {
  // likelihood
  int start_pos;
  row_vector[3] f; // start_vector
  matrix[3,3] P; // transmission probability matrix
  matrix[3,3] E; // emission probability matrix
  start_pos = 1;
  for (i in 1:N) {
    // for participant i
    // setup start vector
    // ( Pr( mic | S(t=0) = i)Pr(S(t=0) = i)) for i = [1,2])
    

    E = generate_E_matrix(state_start[start_pos]);
    f = to_row_vector(p0[i]) * E;

    if (debug == 1) {
      print("participant: ", i);
      //print("mu: ", mu, ". sigma: ", sigma);
      //print("beta parameters: ", beta_parameters);
      //print("p0: ", p0[i]);
      print("p0: ", p0[i]);
      print("E matrix: ", E);
      print("f vector: ", f);
    }

    for (j in start_pos:(n_obs[i] + start_pos - 1)) {
      // generate matrix P
      P = generate_P_matrix(
        t_end[j] - t_start[j],
        q0, q1, q2,
        beta_01, beta_02,
        beta_10, beta_12,
        beta_20, beta_21,
        covs[j],
        debug);
      E = generate_E_matrix(
        state_end[j]);

      f = f * P * E;

      if (debug == 1) {
        print("current observation index: ", j);
        print("start t: ", t_start[j], ". end t: ", t_end[j], ". covs: ", covs[j]);
        print("P matrix:", P);
        print("E matrix:", E);
        print("end f vector: ", f);
      }
    }
    if (debug == 1) {
      print("log-likelihood: ", log(sum(f)));
    }
    target += log(sum(f));
    start_pos = start_pos + n_obs[i];
  }
  // priors
  q0 ~ normal(0.2, 0.2);
  q1 ~ normal(0.2, 0.2);
  q2 ~ normal(0.2, 0.2);
  beta_01 ~ normal(0, 2);
  beta_02 ~ normal(0, 2);
  beta_10 ~ normal(0, 2);
  beta_12 ~ normal(0, 2);
  beta_20 ~ normal(0, 2);
  beta_21 ~ normal(0, 2);
  for (k in 1:N) {
    p0[k] ~ dirichlet(rep_vector(0.01, 3));
  }
// mu[1] ~ normal(0, 1);
// mu[2] ~ normal(-4.2,1);
// sigma ~ normal(0,0.5);
// p0 ~ uniform(0,1);
  //p0 ~ beta(0.5,0.5);
  
  //p0 ~ beta(beta_parameters[1], beta_parameters[2]);
  //beta_parameters ~ normal(0,3);
}

generated quantities {
   // ... declarations ... statements ...
}
