

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
   // initial state
   vector<lower=0, upper=1>[N] p_0; // (pr(S[t=0]) = j, j = [0,1], for participant i
   vector<lower=0>[2] beta_params; // beta distribution parameters
   // state transition
   vector<lower=0>[2] q; // state transition instantaneous probs (q01, q10)
   vector[1] beta_01; // covariates on state transitions q00, q01
   vector[1] beta_10; //
   //mic distributions
   vector[2] mu; // mean of MIC distributions for state (0,1)
   vector<lower=0>[2] sigma; // sd of MIC distributions for state (0,1)
}

transformed parameters {
}

model {
  // likelihood
  int start_pos;
  row_vector[2] f; // start_vector
  matrix[2,2] T;
  start_pos = 1;
  for (i in 1:N) {
    // for participant i
    // setup start vector
    // ( Pr( mic | S(t=0) = i)Pr(S(t=0) = i)) for i = [1,2])
    
    f[1] = p_0[i];
    f[2] = 1 - p_0[i];
    if (mic_start[start_pos] > 0) {
      f[1] = f[1] * exp(normal_lpdf(log(mic_start[start_pos]) | mu[1], sigma[1]));
      f[2] = f[2] * exp(normal_lpdf(log(mic_start[start_pos]) | mu[2], sigma[2]));
    }
    if (debug == 1) {
    print("participant: ", i);
    print("mu: ", mu);
    print("sigma: ", sigma);
    print("q" , q);
    print("mic_start: ", mic_start[start_pos]);
    print("p_0: ", p_0[i]);
    if (mic_start[start_pos] > 0) {
    print("Pr(mic | S = 0) :", exp(normal_lpdf(log(mic_start[start_pos]) | mu[1], sigma[1])));
    print("Pr(mic | S = 1) :", exp(normal_lpdf(log(mic_start[start_pos]) | mu[2], sigma[2])));
    }
    print("f vector: ", f);
    }
    for (j in start_pos:(n_obs[i] + start_pos - 1)) {
      // generate matrix Trs
      T = generate_P_matrix(
        t_end[j] - t_start[j],
        q,
        beta_01,
        beta_10,
        covs[j]);
      if (mic_end[j] > 0) {
        T[,1] = T[,1] * exp(normal_lpdf(log(mic_end[j]) | mu[1], sigma[1]));
        T[,2] = T[,2] * exp(normal_lpdf(log(mic_end[j]) | mu[2], sigma[2]));
      }
      if (debug == 1) {
      print("current observation index: ", j);
      print("start t: ", t_start[j], ". end t: ", t_end[j], ". covs: ", covs[j]);
      print("mic start: ", mic_start[j], "mic end: ", mic_end[j]);
      print(" t matrix:", T);
      }
      f = f * T;
    }
    target += log(sum(f));
    start_pos = start_pos + n_obs[i];
  }
  // p_0
  p_0 ~ beta(beta_params[1], beta_params[2]);
  beta_params ~ normal(0,3);
  // priors
  q ~ student_t(3, 0, 1);
  beta_01 ~ student_t(3, 0, 3);
  beta_10 ~ student_t(3, 0, 3);
  mu ~ normal(0,4);
  sigma ~ normal(0,1);
}

generated quantities {
   // ... declarations ... statements ...
}
