//
// This Stan Ecological Inference Program
// Constrains to row and column margins using: sequential sampling
// with various models of row to column rate

functions{
  #include include/allocationfuns.stan
  #include include/realpdf.stan
  #include include/lkjonionfun.stan

  vector simplex_constrain_softmax_lp(vector v) {
     int K = size(v) + 1;
     vector[K] v0 = append_row(0, v);
     return softmax(v0);
  }

}
data{
 int<lower=0> n_areas;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 matrix<lower=0>[n_areas, R] row_margins; // the row margins in each area
 matrix<lower=0>[n_areas, C] col_margins; // the column margins in each area
 int<lower=0, upper=1> structural_zeros[n_areas, R, C];  // an array indicating any structural zeros in the data (may include whole rows, whole columns and/or individual cells)
 int<lower=0, upper=2> lflag_dist; // flag indicating whether to use poisson (0), multinomial (1) or negative binomial (2) paramertization
 int<lower=0, upper=3> lflag_area_re; // flag indicating whether the area mean simplex is uniform (0) or varies with area random effects which are normally distributed (1) or varies with area random effects which are multinormally distributed (non centred paramaterisation) (2) or varies with area random effects which are multinormally distributed (non centred LKJ Onion paramaterisation)
 int<lower  =0, upper=2> lflag_vary_sd; // flag indicating whether variance of area_cell parameters is: (0) shared across cells,  (1) varies by cell,  or (2) has a hierarchical model structure
 int<lower = 0, upper = 1> lflag_llmod_const; // flag indicating whether log-linear model should be (0) unconstrained (1) constrained to generate correct area totals
 int<lower = 0, upper = 1> lflag_llmod_omit_jr; // flag indicating whether log-linear model should omit area * row interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jc; // flag indicating whether log-linear model should omit area * col interaction
 int<lower = 0, upper = 1> lflag_llmod_omit_jrc; // flag indicating whether log-linear model should omit area * row * column interaction
  int<lower = 0, upper =1> lflag_llmod_centred; // flag indicating whether to use (1) centred or (0) decentred parameterization
 int<lower = 0, upper =2> lflag_llmod_structure; // flag indicating whether log-linear model should be (0) saturated (1) omit area cell effect  (areaxrowxcolumn interaction) (2) omit area*column interaction (but include areaxrowXcolumn effect) [replaced by other flags]
 real<lower=0> prior_lkj; // lkj param
 real<lower=0> prior_mu_re_scale; // prior for scale of mu_re (mean row effect)
 real<lower=0> prior_mu_ce_scale; // prior for scale of mu_ce (mean column effect)
 real<lower=0> prior_sigma_c_scale; //prior for scale of sigma_c (or sigma_c_sigma if lflag_vary_sd == 3)
 real<lower=0> prior_sigma_c_mu_scale; //prior for scale of sigma_c_mu (only if lflag_vary_sd == 3)
 real<lower=0> prior_sigma_ce_scale; //prior of scale for sigma_ce
 real<lower=0> prior_sigma_re_scale; //prior of scale for sigma_re
 real<lower=0> prior_cell_effect_scale; //prior of scale for average cell effects
}
transformed data{
  int K;
  int K_t;
  int K_no_rm;
  int K_c;
  int K_ame;
  int has_area_cell_effects;
  int free_R[n_areas];
  int free_C[n_areas];
  int non0_rm;
  int non0_cm;
  int n_poss_cells=0;
  int n_structural_zeros=0;
  int structural_zero_rows[n_areas, R];
  int n_poss_rows=0;
  int structural_zero_cols[n_areas, C];
  int n_poss_cols=0;
  int has_phi;
  int has_area_re;
  int has_free_area;
  int has_area_col_effects;
  int has_area_row_effects;
  int has_L;
  int has_onion;
  int has_L_ame;
  // int mu_re_ce_in_cell_effects  = 1;
  real param_map[n_areas, R - 1, C - 1];
  matrix[n_areas, R - 1] row_margins_lr;
  real<lower=0> prior_phi_scale;
  int n_margin_sigmas;
  int n_cell_sigmas;
  int n_table_sigmas;

  prior_phi_scale =3;

  if(lflag_dist==2){
    has_phi = 1;
  } else{
    has_phi = 0;
  }

  if(lflag_area_re == 0){
    has_area_re = 0;
  } else {
    has_area_re = 1;
  }

  if(lflag_area_re == 2){
    has_L = 1;
  } else {
    has_L = 0;
  }

  if(lflag_area_re == 3){
    has_onion = 1;
  } else {
    has_onion = 0;
  }

  if(lflag_llmod_const == 0){
    has_free_area = 1;
  } else {
    has_free_area = 0;
  }

  if(lflag_llmod_omit_jrc ==1){
    has_area_cell_effects = 0;
  } else {
    has_area_cell_effects = 1;
  }

  if(lflag_llmod_omit_jc == 1){
    has_area_col_effects = 0;
  } else {
    has_area_col_effects = 1;
  }
  if(lflag_llmod_omit_jr == 1){
    has_area_row_effects = 0;
  } else {
    has_area_row_effects = 1;
  }



  K_t = (R - 1)* (C - 1);
  K_ame = has_area_row_effects * (R - 1) + has_area_col_effects * (C - 1); // number of area margin effects to estimate (area row effects + area col effects)

  K_c = C * (R - 1);
  K = R * (C - 1);
  K_no_rm = R * (C - 1);

  // if(lflag_inc_rm == 1){
  //   K = (R * C) - 1;
  // } else {
  //   K = R * (C - 1);
  // }
  // K_no_rm = R * (C - 1);

  array[2] vector[K -1] shapes = create_shapes(K_t, prior_lkj);

// following code to distinguish structural zeros from sampling zeros in the log-linear model - identifying rows and columns, and number of structural zeros in the data

  for(j in 1:n_areas){
    for(r in 1:R){
      structural_zero_rows[j, r] = sum(structural_zeros[j, r, 1:C])==C ? 1 :0;
      n_poss_rows += sum(structural_zeros[j, r, 1:C])==C ? 0 :1;
      n_structural_zeros += sum(structural_zeros[j, r, 1:C]);
    }
    for(c in 1:C){
      structural_zero_cols[j, c] = sum(structural_zeros[j, 1:R, c])==R ? 1 : 0;
      n_poss_cols += sum(structural_zeros[j, 1:R, c])==R ? 0 :1;
    }
  }
  n_poss_cells = n_areas*R*C - n_structural_zeros;


  // to deal with zeros in the sequential sampling: calculate the number of free parameters (zero row and columns do not need a parameter to allocated cell value of 0)

  int n_param = 0;
  int param_count_from[n_areas];

  for(j in 1:n_areas){
    free_R[j] = 0;
    free_C[j] = 0;
    for(r in 1:R){
      if(row_margins[j, r]>0){
        free_R[j] += 1;
      }
    }
    for(c in 1:C){
      if(col_margins[j, c] > 0){
        free_C[j] += 1;
      }
    }
    param_map[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j] - 1)){
      for (c in 1:(free_C[j] - 1)){
        param_map[j, r, c] = n_param + ((r - 1) * (free_C[j] -1)) + c;
       }
     }
    param_count_from[j] = n_param;
    n_param += max(0, (free_R[j] - 1) * (free_C[j] - 1));
  }
  non0_rm = sum(free_R);
  non0_cm = sum(free_C);

  for (j in 1:n_areas){
    for(r in 1:(R - 1)){
      row_margins_lr[j, r] = log((row_margins[j, r]+ .5)/(sum(row_margins[r, (r + 1):R]) + .5));
    }
  }

  has_L_ame = 0;

  n_margin_sigmas = K_ame;
  n_cell_sigmas = has_area_cell_effects * ((lflag_vary_sd ==0) ? 1 : (R - 1) * (C - 1));
  n_table_sigmas = n_margin_sigmas + n_cell_sigmas;




}
parameters{
 real lambda_unpadded[n_param]; // sequential cell weights
 // vector[lflag_mod_cols * K_c] mu_c;  // logit probability row mean to construct theta

 // matrix[has_area_re*n_areas, has_area_re * R * (C - 1)] alpha; // area variation from mu (logit probability row mean)
 // array[has_area_re*n_areas] vector[has_area_re*K] alpha;    // logit column probabilities for area j and rows
 array[has_area_cell_effects*(n_areas -1)] vector[(R - 1) * (C - 1)] area_cell_effect_raw;    // logit column probabilities for area j and cols
 real area_effect_raw[has_free_area*(n_areas)];
 real area_effect_mu[has_free_area];
 // matrix[n_areas, (R - 1)] area_row_effect_raw;
 // matrix[has_area_col_effects*n_areas, has_area_col_effects*(C - 1)] area_col_effect_raw;
 // array[n_areas] vector[(R - 1) + has_area_col_effects* (C - 1)] ame_alpha;
 array[n_areas] row_vector[K_ame] area_margin_effects_raw;
 // matrix[(R - 1), (C - 1)] cell_effect_raw;
 vector[(R -1) * (C - 1)] mu_area_cell_effect_raw;
 vector[(R -1) * (C - 1)] cell_effect_raw;
 // real<lower=0> cell_effect_sigma;

 // matrix[R, C] cell_effect;
 vector[(C - 1)] mu_ce_raw;
 vector[(R - 1)] mu_re_raw;
 // vector[C - 1] col_effect_raw;
 // vector[R - 1] row_effect_raw;
 // vector<lower=0>[has_area_col_effects*(C - 1)] sigma_ce;
 // vector<lower=0>[R - 1] sigma_re;
 real<lower=0> sigma_area_effect[has_free_area];

 // vector<lower=0>[(R - 1) * (C - 1)] sigma_c; //scale of area col variation from mu_c
 real<lower=0> all_table_sigmas_raw[n_table_sigmas];

  // real<lower=0> sigma_c_raw_raw[(lflag_vary_sd ==0) ? 1 : (R - 1) * (C - 1)]; //scale of area col variation from mu_c
 real<lower=0> sigma_c_sigma[has_area_cell_effects * ((lflag_vary_sd ==2) ? 1 : 0)]; //hierarchical standard deviatation on standard deviations
 vector[has_area_cell_effects * ((lflag_vary_sd ==2) ? 1 : 0)] sigma_c_mu; //hierarchical mean on standard deviations
 // cholesky_factor_corr[has_L * K] L_a; // for modelling correlation matrix
 // cholesky_factor_corr[has_L_ame*K_ame] L_ame_raw; // for modelling correlation matrix
 // array[n_areas, R] simplex[C] theta_jr;
 // array[R] simplex[C] d_theta_r;
 // row_vector[has_onion * (choose(K, 2) - 1)] l; // do NOT init with 0 for all elements
 // vector<lower = 0, upper = 1>[has_onion * (K - 1)] R2; // first element is not really a R^2 but is on (0,1)
 // matrix[lflag_predictors_rm * K_no_rm, lflag_predictors_rm * (R - 1)] betas_rm;
 real mu_cell_effects_raw;
 real<lower=0> scale_cell_effects;
 real mu_col_effects_raw;
 real<lower=0> scale_col_effects;
 real mu_row_effects_raw;
 real<lower=0> scale_row_effects;
 real mu_effects;
 real<lower=0> scale_effects;
 real<lower=0, upper=1> phi[has_phi];
 real mu_log_phi[has_phi];
 real<lower=0> log_phi_scale[has_phi];
 // real overall_mu;
 real<lower=0, upper = .001> hinge_delta_floor;
 real<lower=0, upper = .001> hinge_delta_min;
}
transformed parameters{
  real lambda[n_areas, R - 1, C -1]; // sequential cell weights
  real<lower=0> cell_values[n_areas, R, C];
  // array[n_areas, R] simplex[C] theta_area; // prob row vector for each area
  // array[n_areas, C] simplex[R] theta_c_area; // prob col vector for each area
  array[n_areas] matrix[R, C] area_effect_components;
  real log_e_cell_value[n_areas, R, C];
  // matrix[n_areas, (R - 1) + has_area_col_effects* (C - 1)] area_margin_effects_raw;

  vector[R] mu_re;
  vector[C] mu_ce;
  vector[K_ame] mu_ame;
  matrix[n_areas, R] area_row_effect;
  matrix[n_areas, C] area_col_effect;
  vector[n_areas] area_effect;
  matrix[n_areas, (R - 1)] area_row_effect_raw;
  matrix[n_areas, (C - 1)] area_col_effect_raw;
  vector[(R -1) * (C - 1)] mu_area_cell_effect;
  array[n_areas] matrix[R, C] area_cell_effect;

  real<lower=0> sigma_c[has_area_cell_effects * (R - 1) * (C - 1)]; //scale of area col variation from mu_c
  real mu_cell_effects;
  real mu_row_effects;
  real mu_col_effects;
  real<lower=0> sigma_c_raw[n_cell_sigmas]; //scale of area col variation from mu_c
  vector<lower=0>[n_margin_sigmas] ame_sigma;
  // real hinge_delta_floor = 1/(100 + inv_hinge_delta_floor);
  // real hinge_delta_min = 1/(100 + inv_hinge_delta_min);
  // // matrix[has_L_ame * K_ame, has_L_ame * K_ame] L_ame;

  // matrix[max(has_onion, has_L) * K, max(has_onion, has_L) * K] L;

  // if(has_onion == 1){
  //   L = lkj_onion(K, l, R2, shapes); // cholesky_factor corr matrix
  // }
  // if(has_L == 1){
  //   L = L_a; // cholesky_factor corr matrix
  // }

  // if(has_L_ame == 1){
  //   L_ame = L_ame_raw;
  // } else{
  //   L_ame = diag_matrix(rep_vector(1, K_ame));
  // }



  for (j in 1:n_areas){
    lambda[j] = rep_array(0, R - 1, C - 1);
    for (r in 1:(free_R[j]-1)){
      for (c in 1:(free_C[j] - 1)){
        lambda[j, r, c] = lambda_unpadded[param_count_from[j] + ((r - 1) * (free_C[j] - 1)) + c];
     }
   }
 }


  cell_values = ss_assign_cvals_wzeros_hinge_lp(n_areas, R, C, row_margins, col_margins, lambda, hinge_delta_floor, hinge_delta_min);

  sigma_c_raw = all_table_sigmas_raw[1:n_cell_sigmas];
  ame_sigma = to_vector(all_table_sigmas_raw[n_cell_sigmas + 1:n_cell_sigmas + n_margin_sigmas]);

  if(lflag_llmod_centred == 1){
    mu_cell_effects = mu_cell_effects_raw;
    mu_col_effects = mu_col_effects_raw;
    mu_row_effects = mu_row_effects_raw;
  } else if(lflag_llmod_centred == 0){
    mu_cell_effects = mu_effects + scale_effects*mu_cell_effects_raw;
    mu_col_effects = mu_effects + scale_effects*mu_col_effects_raw;
    mu_row_effects = mu_effects + scale_effects*mu_row_effects_raw;
  }

  mu_re = rep_vector(0, R);
  mu_ce = rep_vector(0, C);
  area_effect = rep_vector(0, n_areas);

    if(lflag_llmod_centred==1){
      mu_re[1:R -1] = mu_re_raw;
    } else if (lflag_llmod_centred==0){
      mu_re[1:R - 1] = mu_row_effects + scale_row_effects * mu_re_raw;
    }
    if(has_area_row_effects == 1){
        mu_ame[1:R - 1] = mu_re[1:R - 1];
    }



    if(lflag_llmod_centred == 0){
      mu_ce[1: C - 1] = mu_col_effects + scale_col_effects*mu_ce_raw;
    } else if(lflag_llmod_centred == 1){
      mu_ce[1: C - 1] = mu_ce_raw;
    }
    if(has_area_col_effects == 1){
      mu_ame[has_area_row_effects*(R - 1) + 1:has_area_row_effects*(R - 1) + 1 + C - 2] = mu_ce[1:C - 1];
    }

    if(lflag_llmod_centred == 0){
      mu_area_cell_effect = mu_cell_effects + scale_cell_effects*mu_area_cell_effect_raw;
    } else if(lflag_llmod_centred ==1){
      mu_area_cell_effect = mu_area_cell_effect_raw;
    }


  // for(r in 1:R){
  //   for(c in 1:C){
  //     if(r == R && c ==C){
  //       cell_effect[r, c] = 0;
  //     } else{
  //       cell_effect[r, c] = cell_effect_raw[(r - 1)*C + c];
  //     }
  //
  //   }
  // }


  for(j in 1:n_areas){

    area_cell_effect[j] = rep_matrix(0.0, R, C);
    for(r in 1:(R - 1)){
      if(has_area_cell_effects == 1 && j < n_areas){
        area_cell_effect[j, r, 1:(C - 1)] = to_row_vector(area_cell_effect_raw[j, ((r - 1)*(C - 1) + 1):r*(C - 1)]);
      } else if(has_area_cell_effects == 0){
        area_cell_effect[j, r, 1:(C - 1)] = to_row_vector(mu_area_cell_effect[((r - 1)*(C - 1) + 1):r*(C - 1)]);
      }

    // } else if(has_area_cell_effects==0){
    //   for(r in 1:(R - 1)){
    //     for(c in 1:(C  - 1)){
    //       area_cell_effect[j, r, c] = mu_area_cell_effect[((C - 1) * (r - 1)) + c ];
    //     }
    //   }
    }

    // area_margin_effects_raw[j] = to_row_vector(mu_ame + ame_sigma .*(ame_alpha[j]));

    area_row_effect_raw[j] = rep_row_vector(0, R - 1);
    area_row_effect[j] = rep_row_vector(0, R);

   if(has_area_row_effects == 1){
      area_row_effect_raw[j] = area_margin_effects_raw[j, 1:(R - 1)];
      area_row_effect[j, 1:(R - 1)] = area_row_effect_raw[j];
    } else if(has_area_row_effects==0){
      area_row_effect_raw[j] = to_row_vector(mu_re[1:(R - 1)]);
      area_row_effect[j, 1:(R - 1)] = area_row_effect_raw[j];
    }


    area_col_effect_raw[j] = rep_row_vector(0, C - 1);
    area_col_effect[j] = rep_row_vector(0, C);

    if(has_area_col_effects == 1 && j < n_areas){
      area_col_effect_raw[j, 1:(C - 1)] = area_margin_effects_raw[j, (R - 1 + 1): (R - 1 + C - 1)];
      area_col_effect[j, 1:(C - 1)] = area_col_effect_raw[j];
    } else if(has_area_col_effects==0){
      area_col_effect_raw[j] = to_row_vector(mu_ce[1:(C - 1)]);
      area_col_effect[j, 1:(R - 1)] = area_col_effect_raw[j];
    }


    for(r in 1:R){
      for(c in 1:C){
        area_effect_components[j, r, c] = area_row_effect[j, r] + area_col_effect[j, c] + area_cell_effect[j, r, c];
      }
    }
  if(has_free_area == 1){
    // if(j<n_areas){
      if(lflag_llmod_centred==1){
        area_effect[j] = area_effect_raw[j];
      } else if (lflag_llmod_centred == 0){
        area_effect[j] = area_effect_mu[1] + sigma_area_effect[1]*area_effect_raw[j];
      }
    // }
  } else {
    area_effect[j] = log(sum(row_margins[j])) - log_sum_exp(to_vector(area_effect_components[j]));
  }

    for(r in 1:R){
        for(c in 1:C){
          log_e_cell_value[j, r, c] = area_effect[j] + area_row_effect[j, r] + area_col_effect[j, c] + area_cell_effect[j, r, c];
        }
      }
    }

    if(has_area_cell_effects==1){
      if(lflag_vary_sd == 0){
        sigma_c = rep_array(sigma_c_raw[1], (R - 1)*(C - 1));
      } else {
          sigma_c = sigma_c_raw;
      }
    }


  #include include/generateratesandsummaries.stan


}
model{
  // matrix[non0_rm, C] obs_prob;
  // vector[n_poss_rows] e_rm;
  // row_vector[n_poss_rows] poss_rows;
  // vector[n_poss_cols] e_cm;
  // row_vector[n_poss_cols] poss_cols;
  vector[n_poss_cells] e_cell;
  vector[has_phi*n_poss_cells] alpha_cell;
  row_vector[n_poss_cells] cell_values_row_vector;
  // matrix[phi_length*non0_rm, phi_length*C] phi_matrix;
  // matrix[non0_cm, R] cell_values_matrix_c;
  // int counter_c = 1;
  // matrix[non0_cm, R] obs_prob_c;
  array[n_areas] vector[C - 1] mu_col_effect_raw;
  array[n_areas] vector[R - 1] mu_row_effect_raw;

  // vector[R * C - 1] mu_cell_effect;

  // if(has_L_ame==1){
  //    L_ame_raw ~ lkj_corr_cholesky(prior_lkj);
  // }


  // for(r in 1:R){
  //   d_theta_r[r] ~ dirichlet(rep_vector(1, C));
  // }


//
//   if(has_L == 1){
//     L ~ lkj_corr_cholesky(prior_lkj); // implies L*L'~ lkj_corr(prior_lkj);
//   }
//   if(has_onion == 1){
//     l ~ std_normal();
//     R2 ~ beta(shapes[1], shapes[2]);
//   }

  // int counter_r=0;
  // int counter_c=0;
  int counter_cell=0;


  for (j in 1:n_areas){
    // for (r in 1:R){
    //   if(structural_zero_rows[j,r]==0){
    //     counter_r += 1;
    //     e_rm[counter_r] = sum(exp(log_e_cell_value[j, r, 1:C]));
    //     poss_rows[counter_r] = row_margins[j, r];
    //   }
   // }

     for (c in 1:C){
       // if(structural_zero_cols[j,c]==0){
       //   counter_c += 1;
       //   e_cm[counter_c] = sum(exp(log_e_cell_value[j, 1:R, c]));
       //   poss_cols[counter_c] = col_margins[j,c];
       // }
       for(r in 1:R){
         if(structural_zeros[j,r,c]==0){
           counter_cell += 1;
           cell_values_row_vector[counter_cell] = cell_values[j,r,c];
           if(lflag_dist == 0){
             e_cell[counter_cell] = exp(log_e_cell_value[j, r, c]);
           } else if(lflag_dist == 2){
             alpha_cell[counter_cell] = phi[1];
             e_cell[counter_cell] = exp(log_e_cell_value[j, r, c])/(phi[1] + exp(log_e_cell_value[j, r, c]));
           }
         }
       }
       // if(col_margins[j, c] > 0){
       //   for(r in 1:R){
       //     cell_values_matrix_c[counter_c, r] = cell_values[j, r, c];
       //     // obs_prob_c[counter_c, r] = theta_c_area[j, c, r] * col_margins[j, c];
       //     obs_prob_c[counter_c, r] = exp(log_e_cell_value[j, r, c]);
       //   }
       //   counter_c += 1;
        // }
      }
    }
  // if(lflag_dist == 1){
  //   target += realmultinom_lpdf(cell_values_matrix | obs_prob);
  // } else if(lflag_dist == 0) {
  //   target += realpoisson_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob));
  // } else if(lflag_dist == 2) {
  //   phi ~ cauchy(0, 3);
  //   target += realnegbinom3_lpdf(to_row_vector(cell_values_matrix) | to_vector(obs_prob), to_vector(phi_matrix));
  // }
  // if(lflag_inc_rm){
  //     target += realmultinom_lpdf(row_margins| rm_prob);
  // }


    // target +=realpoisson_lpdf(to_row_vector(cell_values_matrix_c)| to_vector(obs_prob_c));
    if(lflag_dist == 0){
      target +=realpoisson_lpdf(cell_values_row_vector| e_cell);
    } else if(lflag_dist==1){
      phi ~ normal(0, prior_phi_scale);
      mu_log_phi ~ normal(0, 1);
      log_phi_scale ~ normal(0, 1);
      target +=realnegbinom3_lpdf(cell_values_row_vector | alpha_cell, e_cell);
    }

    // target +=realpoisson_lpdf(poss_rows| e_rm);
    // target +=realpoisson_lpdf(poss_cols| e_cm);

    // for(r in 1:R){
    //   for(c in 1:C){
    //     e_overall_cv[C * (r - 1) + c]=sum(exp(log_e_cell_value[1:n_areas, r, c]));
    //   }
    // }
    // vector[R * C] e_overall_cv;
    // target +=realpoisson_lpdf(to_row_vector(overall_cell_values')| e_overall_cv);

    // print(to_row_vector(overall_cell_values'));
    // print(e_overall_cv);


    // sigma_c ~ normal(0, 5);
    if(has_area_cell_effects==1){
      if(lflag_vary_sd == 2){
        sigma_c_mu ~ normal(0, prior_sigma_c_mu_scale);
        sigma_c_sigma ~ normal(0, prior_sigma_c_scale);
        for(s in 1:K_t){
          sigma_c_raw[s]~lognormal(sigma_c_mu, sigma_c_sigma);
        }
      } else {
        sigma_c_raw ~ normal(0, prior_sigma_c_scale);
      }
    }
    // mu_ce_raw~ normal(0, prior_mu_ce_scale);
    // mu_re_raw~ normal(0, prior_mu_re_scale);
    // mu_area_cell_effect_raw ~ normal(0, prior_cell_effect_scale);

    if(lflag_llmod_centred ==0 ){
       mu_ce_raw ~ std_normal();
       mu_re_raw ~ std_normal();
       mu_area_cell_effect_raw ~ std_normal();
    } else if(lflag_llmod_centred==1){
       mu_ce_raw ~ normal(mu_col_effects, scale_col_effects);
       mu_re_raw~ normal(mu_row_effects, scale_row_effects);
       mu_area_cell_effect_raw ~ normal(mu_cell_effects, scale_cell_effects);
    }


    scale_effects ~ normal(0, prior_cell_effect_scale);
    mu_effects ~ normal(0, prior_cell_effect_scale);

    if(lflag_llmod_centred == 1){
      mu_row_effects ~ normal(mu_effects, scale_effects);
      mu_col_effects_raw ~ normal(mu_effects, scale_effects);
      mu_cell_effects ~ normal(mu_effects, scale_effects);
    } else {
      mu_col_effects_raw ~ std_normal();
      mu_row_effects_raw ~ std_normal();
      mu_cell_effects_raw ~ std_normal();
    }

    // mu_cell_effects ~ normal(0, prior_cell_effect_scale);
    scale_cell_effects~normal(0, prior_cell_effect_scale);
    scale_col_effects~normal(0, prior_mu_ce_scale);
    scale_row_effects~normal(0, prior_mu_re_scale);


    // sigma_re~normal(0, prior_sigma_re_scale);
    // to_vector(cell_effect_raw) ~ normal(0, prior_cell_effect_scale);
    // for(r in 1:R){
    //   for(c in 1:C){
    //     if( (r + c) < (R + C)){
    //        mu_cell_effect[(r - 1)*C + c] = log(overall_cell_values[r, c]) - log(overall_cell_values[R, C]);
    //     }
    //   }
    // }
    // // cell_effect_raw ~ normal(mu_cell_effect, cell_effect_sigma);
    // cell_effect_sigma ~ cauchy(0,2);
    // to_vector(cell_effect) ~ normal(0, prior_cell_effect_scale);
    for (j in 1:(n_areas - 1)){
      // ame_alpha[j] ~ std_normal();
      if(has_area_row_effects==1 || has_area_col_effects == 1){
          area_margin_effects_raw[j] ~ normal(mu_ame, ame_sigma);
      }
      if(has_area_cell_effects==1){
          to_vector(area_cell_effect_raw[j]) ~ normal(mu_area_cell_effect_raw, sigma_c);
      }
      // area_row_effect_raw[j] ~ normal(0, sigma_re);
      // if(has_area_col_effects == 1){
      //    area_col_effect_raw[j] ~ normal(0, sigma_ce);
      //    sigma_ce~normal(0, prior_sigma_ce_scale);
      //
      // }

      // area_row_effect_raw[j] ~ normal(0, sigma_re);
      // area_col_effect_raw[j] ~ normal(0, sigma_ce);

      // area_col_effects are to meet column proportion sufficient statistics
      // area_col effect = log(Tc) + log_sum_exp(row_effect + cell_effect [within the reference column]) - log(Tref) - log_wum_exp(row_effect + cell_effect[within the column])


    }
    if(lflag_vary_sd == 2 && has_area_cell_effects == 1){
       ame_sigma~lognormal(sigma_c_mu[1], sigma_c_sigma[1]);
    } else{
       ame_sigma~normal(0, 3);
    }


    if(has_free_area == 1){
      if(lflag_llmod_centred  == 1){
        area_effect_raw ~ normal(area_effect_mu[1], sigma_area_effect[1]);
      } else {
        area_effect_raw ~ std_normal();
      }
      sigma_area_effect~normal(0, 3);
      area_effect_mu~normal(0, 10);
      // overall_mu ~ normal(0, 3);
    }
    hinge_delta_floor ~ normal(0, .001);
    hinge_delta_min ~ normal(0, .001);

}
generated quantities{

}
