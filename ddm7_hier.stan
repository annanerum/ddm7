functions {
  real partial_sum(array[] int idx_slice, int start, int end,
                   vector rt, array[] int acc, array[] int pid, vector condition,
                   vector a, vector t0, vector w, vector sv, vector sw, vector st0,
                   vector nu_intercept, vector nu_c, data real prec) {
    real lp = 0;
    for (n in 1:size(idx_slice)) {
      int i = idx_slice[n];
      int p = pid[i];
      real v = nu_intercept[p] + nu_c[p] * condition[i];
      if (acc[i] == 1) {
        lp += wiener_lpdf(rt[i] | a[p], t0[p], w[p], v, sv[p], sw[p], st0[p], prec);
      } else {
        lp += wiener_lpdf(rt[i] | a[p], t0[p], 1 - w[p], -v, sv[p], sw[p], st0[p], prec);
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> I;
  array[N] int<lower=1, upper=I> pid;
  vector<lower=0>[N] rt;
  array[N] int<lower=-1, upper=1> acc;
  vector[N] condition;
  vector<lower=0>[I] min_rt;
  int<lower=1> grainsize;

  vector[2] prior_beta_nu_intercept;
  vector[2] prior_beta_nu_c;
  vector[2] prior_beta_alpha;
  vector[2] prior_beta_tau;
  vector[2] prior_beta_w;
  vector[2] prior_beta_sv;
  vector[2] prior_beta_sw;                  // logit of the fraction of the max sw
  vector[2] prior_beta_st0;                 // log scale (log seconds)

  vector[2] prior_sigma_nu_intercept;
  vector[2] prior_sigma_nu_c;
  vector[2] prior_sigma_alpha_intercept;
  vector[2] prior_sigma_tau_intercept;
  vector[2] prior_sigma_w;
  vector[2] prior_sigma_sv;
  vector[2] prior_sigma_sw;                 
  vector[2] prior_sigma_st0;                

  
  real guard_a_max;      // max boundary separation
  real guard_sv_max;     // max drift variability
  real guard_v_max;      // max |drift| in any condition
  real guard_w_min;      // w must lie in (guard_w_min, 1 - guard_w_min)
  real guard_st0_max;    // max NDT variability (s)

  real<lower=0> precision_derivatives;   
}

transformed data {
  array[N] int idx;
  for (n in 1:N) idx[n] = n;
}

parameters {
  real beta_nu_intercept;
  real beta_nu_c;
  real beta_alpha;
  real<lower=0.05, upper=0.45> beta_tau;
  real beta_w;
  real beta_sv;
  real beta_sw;                             
  real beta_st0;                            

  real<lower=0> sigma_nu_intercept;
  real<lower=0> sigma_nu_c;
  real<lower=0> sigma_alpha_intercept;
  real<lower=0> sigma_tau_intercept;
  real<lower=0> sigma_w;
  real<lower=0> sigma_sv;
  real<lower=0> sigma_sw;                   
  real<lower=0> sigma_st0;                  

  vector[I] z_nu_intercept;
  vector[I] z_nu_c;
  vector[I] z_alpha_intercept;
  vector<lower=0, upper=min_rt>[I] t0;
  vector[I] z_w;
  vector[I] z_sv;
  vector[I] z_sw;                           
  vector[I] z_st0;                          
}

transformed parameters {
  vector[I] nu_intercept    = beta_nu_intercept + sigma_nu_intercept    * z_nu_intercept;
  vector[I] nu_c            = beta_nu_c         + sigma_nu_c            * z_nu_c;
  vector[I] alpha_intercept = beta_alpha        + sigma_alpha_intercept * z_alpha_intercept;
  vector[I] w_logit         = beta_w            + sigma_w               * z_w;
  vector[I] log_sv          = beta_sv           + sigma_sv              * z_sv;
  vector[I] log_st0         = beta_st0          + sigma_st0             * z_st0;   

  vector[I] a   = exp(alpha_intercept);
  vector[I] w   = inv_logit(w_logit);
  vector[I] sv  = exp(log_sv);
  vector[I] st0 = exp(log_st0);                                                     
  vector[I] sw;                                                                     
  for (i in 1:I)
    sw[i] = 0.99 * 2 * fmin(w[i], 1 - w[i]) * inv_logit(beta_sw + sigma_sw * z_sw[i]);
}

model {
  beta_nu_intercept ~ normal(prior_beta_nu_intercept[1], prior_beta_nu_intercept[2]);
  beta_nu_c         ~ normal(prior_beta_nu_c[1],         prior_beta_nu_c[2]);
  beta_alpha        ~ normal(prior_beta_alpha[1],        prior_beta_alpha[2]);
  beta_tau          ~ normal(prior_beta_tau[1],          prior_beta_tau[2]);
  beta_w            ~ normal(prior_beta_w[1],            prior_beta_w[2]);
  beta_sv           ~ normal(prior_beta_sv[1],           prior_beta_sv[2]);
  beta_sw           ~ normal(prior_beta_sw[1],           prior_beta_sw[2]);
  beta_st0          ~ normal(prior_beta_st0[1],          prior_beta_st0[2]);

  sigma_nu_intercept    ~ normal(prior_sigma_nu_intercept[1],    prior_sigma_nu_intercept[2]);
  sigma_nu_c            ~ normal(prior_sigma_nu_c[1],            prior_sigma_nu_c[2]);
  sigma_alpha_intercept ~ normal(prior_sigma_alpha_intercept[1], prior_sigma_alpha_intercept[2]);
  sigma_tau_intercept   ~ normal(prior_sigma_tau_intercept[1],   prior_sigma_tau_intercept[2]);
  sigma_w               ~ normal(prior_sigma_w[1],               prior_sigma_w[2]);
  sigma_sv              ~ normal(prior_sigma_sv[1],              prior_sigma_sv[2]);
  sigma_sw              ~ normal(prior_sigma_sw[1],              prior_sigma_sw[2]);
  sigma_st0             ~ normal(prior_sigma_st0[1],             prior_sigma_st0[2]);

  z_nu_intercept    ~ std_normal();
  z_nu_c            ~ std_normal();
  z_alpha_intercept ~ std_normal();
  t0 ~ normal(beta_tau, sigma_tau_intercept);
  z_w               ~ std_normal();
  z_sv              ~ std_normal();
  z_sw              ~ std_normal();
  z_st0             ~ std_normal();

  
  if (max(a) > guard_a_max || max(sv) > guard_sv_max || max(st0) > guard_st0_max ||
      min(w) < guard_w_min || max(w) > 1 - guard_w_min ||
      max(abs(nu_intercept) + abs(nu_c)) > guard_v_max) {
    target += negative_infinity();
  } else {
    target += reduce_sum(partial_sum, idx, grainsize,
                         rt, acc, pid, condition,
                         a, t0, w, sv, sw, st0, nu_intercept, nu_c,
                         precision_derivatives);
  }
}

generated quantities {
  real w_pop_median     = inv_logit(beta_w);
  real sv_pop_median    = exp(beta_sv);
  real alpha_pop_median = exp(beta_alpha);
  real st0_pop_median   = exp(beta_st0);
  real sw_pop_mean      = mean(sw);         
  vector[I] mean_ndt    = t0 + st0 / 2;
}
