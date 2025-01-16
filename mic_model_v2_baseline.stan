// random intercept per participant
// Time-varying effect of drug effect  

data {
int<lower=1> n;                 //total number of observations
int n_participants;             // number of participants
array[n] int pid;                     // map observation on to participant id
vector[n] t_e;                  // days since drug exposure; <0 = no exposure
vector[n] y;                    // response variable (log MIC)
vector[n] y_0;      // value at baseline
}


parameters{
// ----------parameters describing MIC distribution --------------------------
real<lower=0> sigma_mic;             // sd of mic distribution
real mu0;                         // mean of mic distruibution
//-----------participant effects ---------------------------------------------
vector[n_participants] gamma;    // intercept for each participant
real<lower=0> sigma_participant;       // sd of participant intercepts
// ----------species and drug effects ----------------------------------------
vector[1] beta;                  // effect on mic of drug exposure
real<lower=0> tau;            // mean lifetime of drug effect
}

transformed parameters{
// transformations to make clear what is doing what
vector[n] ind_effect;            // individual effect for each sample
vector[n] drug_effect;           // drug effect for each sample

for (i in 1:n) {
  if(t_e[i] < 0) {
    drug_effect[i] = 0;
  }
  else {
    drug_effect[i] = beta[1]*exp(-t_e[i] / tau);
  }
  ind_effect[i] = gamma[pid[i]];
}

vector[n] mu; // mean for normal log(MIC) distribution for each sample
mu = mu0 + y_0 + drug_effect + ind_effect; 
}

model{
// likelihood -----------------------------------------------------------------
y ~ normal(mu, sigma_mic);
// normally distributed individual effects ------------------------------------
gamma ~ normal(0, sigma_participant);
//priors - weakly informative -------------------------------------------------
sigma_mic ~ student_t(3, 0, 3);
beta ~ student_t(3, 0, 3);
sigma_participant ~ student_t(3, 0, 3);
tau ~ student_t(3,0,30);
mu0 ~ student_t(3,0,3);
}

generated quantities {
// save predicted y from the model for posterior checking ---------------------
array[n] real y_pred = lognormal_rng(mu, sigma_mic);
}

