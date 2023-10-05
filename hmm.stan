

functions {

  matrix generate_P_matrix(real t, vector q, vector beta_01, vector beta_10, vector covs) {
    matrix[2,2] Q;
    matrix[2,2] P;

    Q[1,2] = q[1] * exp(dot_product(beta_01, covs));
    Q[1,1] = -Q[1,2];
    Q[2,1] = q[2] * exp(dot_product(beta_10, covs));
    Q[2,2] = -Q[2,1];
    P = matrix_exp(Q * t);
    return(P);
  }

  matrix generate_E_matrix(real mic, vector mu, vector sigma) {
    matrix[2,2] E;
    
    if (mic > 0) {
      E[1,1] = exp(normal_lpdf(log(mic) | mu[1], sigma[1]));
    } else {
      E[1,1] = 1; 
    }
    
    E[1,2] = 0;
    E[2,1] = 0;
    
    if (mic > 0) {
      E[2,2] = exp(normal_lpdf(log(mic) | mu[2], sigma[1]));
    } else {
      E[2,2] = 1;
    
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
  // time and mic at start and end of observation-pair
  // negative mic value encodes censored value
  array[n] real t_start; 
  array[n] real t_end;
  array[n] real mic_start;
  array[n] real mic_end;
  // covariate values for segment
  array[n] vector[1] covs;

}

transformed data {
   // ... declarations ... statements ...
}

parameters {
   // state transition
   vector<lower=0>[2] q; // state transition instantaneous probs (q01, q10)
   vector[1] beta_01; // covariates on state transitions q00, q01
   vector[1] beta_10; //
   //mic distributions
   vector[2] mu; // mean of MIC distributions for state (0,1)
   vector<lower=0>[1] sigma; // sd of MIC distributions for state (0,1)
   // initial states
   vector<lower=0, upper=1>[N] p0;
   vector<lower=0>[2] beta_parameters;
}

transformed parameters {
}

model {
  // likelihood
  int start_pos;
  row_vector[2] f; // start_vector
  matrix[2,2] P; // transmission probability matrix
  matrix[2,2] E; // emission probability matrix
  start_pos = 1;
  for (i in 1:N) {
    // for participant i
    // setup start vector
    // ( Pr( mic | S(t=0) = i)Pr(S(t=0) = i)) for i = [1,2])
    f = [p0[i], 1 - p0[i]];
    E = generate_E_matrix(mic_start[start_pos], mu, sigma);
    f = f * E;

    if (debug == 1) {
      print("participant: ", i);
      print("mu: ", mu, ". sigma: ", sigma);
      print("q" , q);
      print("mic_start: ", mic_start[start_pos]);
      print("beta parameters: ", beta_parameters);
      print("p0: ", p0[i]);
      print("f vector: ", f);
      print("E matrix: ", E);
    }

    for (j in start_pos:(n_obs[i] + start_pos - 1)) {
      // generate matrix P
      P = generate_P_matrix(
        t_end[j] - t_start[j],
        q,
        beta_01,
        beta_10,
        covs[j]);
      E = generate_E_matrix(
        mic_end[j],
        mu,
        sigma);


      f = f * P * E;

      if (debug == 1) {
        print("current observation index: ", j);
        print("start t: ", t_start[j], ". end t: ", t_end[j], ". covs: ", covs[j]);
        print("mic start: ", mic_start[j], ". mic end: ", mic_end[j]);
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
  q ~ normal(0, 0.5);
  beta_01 ~ student_t(3, 0, 2);
  beta_10 ~ student_t(3, 0, 2);
  mu ~ normal(0,4);
  sigma ~ normal(0,4);
  p0 ~ beta(beta_parameters[1], beta_parameters[2]);
  beta_parameters ~ normal(0,3);
}

generated quantities {
   // ... declarations ... statements ...
}
