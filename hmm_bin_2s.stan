

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

  matrix generate_E_matrix(int state) {
    matrix[2,2] E;
    E= rep_matrix(0, 2, 2);
    if (state < 0) {
    E[1, 1] = 1;
    E[2, 2] = 1;
    } else {
      if (state == 0) {
        E[1,1] = 1;
      } else { 
        E[2,2] = 1;
      }
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
   vector<lower=0>[2] q; // state transition instantaneous probs (q01, q10)
   vector[1] beta_01; // covariates on state transitions q00, q01
   vector[1] beta_10; //
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
    if (state_start[start_pos] == 0) {
      f = [1, 0];
    } else {
      f = [0, 1];
    }
    // Only needed for hmm
    //E = generate_E_matrix(mic_start[start_pos], breakpoint);
    //f = f * E;

    if (debug == 1) {
      print("participant: ", i);
      print("q" , q);
      print("state_start: ", state_start[start_pos]);
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
        state_end[j]);


      f = f * P * E;

      if (debug == 1) {
        print("current observation index: ", j);
        print("start t: ", t_start[j], ". end t: ", t_end[j], ". covs: ", covs[j]);
        print("state start: ", state_start[j], ". state end: ", state_end[j]);
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
  q ~ normal(0.1, 0.5);
  beta_01 ~ student_t(3, 0, 3);
  beta_10 ~ student_t(3, 0, 3);
}

generated quantities {
   // ... declarations ... statements ...
}
