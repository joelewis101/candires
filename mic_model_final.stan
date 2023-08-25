// normally distributed log(mic) model
// random intercept per participant
// random intercept for species
// random slope by species for drug effect
// Time-varying effect of drug effect  

data {
int<lower=1> n;                 //total number of observations
int n_participants;             // number of participants
int n_species;                  // number of Candida species
int pid[n];                     // map observation on to participant id
int species_id[n];              // map observation onto species
vector[n] t_e;                  // days since drug exposure; <0 = no exposure
vector[n] y;                    // response variable (log MIC)
}


parameters{
// ----------parameters describing MIC distribution --------------------------
real<lower=0> sigma;             // sd of mic distribution
//-----------participant effects ---------------------------------------------
vector[n_participants] gamma;    // intercept for each participant
real<lower=0> gamma_sigma;       // sd of participant intercepts
// ----------species and drug effects ----------------------------------------
array[n_species] vector[2] beta; // a length 2 vector for each species
                                 // encoding the effect of that species on MIC
                                 // distribution (beta1 = intercept)
                                 // and the effect of drug on that species MIC
                                 // distribution (beta2 = slope)
vector<lower=0>[2] beta_sigma;   // sd of beta1 and beta2 distribtions
corr_matrix[2] Omega;            // correlation matrix for intercept and slope
real<lower=0> tau;            // mean lifetime of drug effect; constant
                                 // across all species
}

transformed parameters{
// transformations to make clear what is doing what
vector[n] ind_effect;            // individual effect for each sample
vector[n] drug_effect;           // drug effect for each sample
vector[n] spec_effect;           // species effect for each sample

for (i in 1:n) {
  if(t_e[i] < 0) {
    drug_effect[i] = 0;
  }
  else {
    drug_effect[i] = exp(-t_e[i] / tau);
  }
  ind_effect[i] = gamma[pid[i]];
  spec_effect[i] = beta[species_id[i], 1] + beta[species_id[i], 2]*drug_effect[i];
}

vector[n] mu; // mean for normal log(MIC) distribution for each sample
mu = spec_effect + ind_effect; 
}

model{
// likelihood -----------------------------------------------------------------
y ~ normal(mu, sigma);
// multivariate normally distributed species, and species-drug effects --------
beta ~ multi_normal(rep_vector(0,2), quad_form_diag(Omega, beta_sigma));
// normally distributed individual effects ------------------------------------
gamma ~ normal(0, gamma_sigma);
//priors - weakly informative -------------------------------------------------
sigma ~ student_t(3, 0, 6);
beta_sigma ~ student_t(3, 0, 6);
Omega ~ lkj_corr(2);
gamma_sigma ~ student_t(3, 0, 6);
tau ~ student_t(3,0,200);
}

generated quantities {
// save predicted y from the model for posterior checking ---------------------
real y_pred[n] = normal_rng(mu, sigma);
}

