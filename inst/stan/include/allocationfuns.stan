
      real robust_hinge_min(vector x, real delta){
        return(-1*delta*log_sum_exp(-1 * x/delta));
      }

      real robust_hinge_floor_zero(real x, real delta) {
        return(delta*log1p_exp(x/delta));
      }

real[,,] ss_assign_lor_wzeros_hinge_newparam_lp(
    int n_areas, int R, int C,
    matrix row_margins, matrix col_margins,
    real[,,] lambda,
    array[] real lambda_zero,
    array[,,] int zero_cell_map,
    array[,,] int structural_zeros,
    real delta_floor, real delta_min) {

     // Declare variables
     real log_cv_out[n_areas, R, C]; // The new return array
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R[n_areas];
     int free_C[n_areas];
     real this_inv_logit;
     real log_det_J = 0.0;

     lower_pos[1] = 0.0;

     // =========================================================================
     // 1. SEQUENTIAL SAMPLING ENGINE
     // =========================================================================
     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);
       free_R[j] = 0;
       free_C[j] = 0;

       for (r in 1:R){
         if(row_margins[j, r] > 0){
           free_R[j] += 1;
           slack_row_raw[free_R[j]] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c] > 0){
           free_C[j] += 1;
           slack_col_raw[free_C[j]] = col_margins[j, c];
         }
       }
       row_vector[free_R[j]] slack_row = slack_row_raw[1:free_R[j]];
       row_vector[free_C[j]] slack_col = slack_col_raw[1:free_C[j]];
       rt = sum(slack_row);
       matrix[free_R[j], free_C[j]] tmp_cell_value;

       for (r in 1:(free_R[j] - 1)){
         for (c in 1:(free_C[j] - 1)){
           lower_pos[2] = slack_row[r] - sum(tail(slack_col, free_C[j] - c));
           lower_bound = robust_hinge_floor_zero(lower_pos[2], delta_floor);
           upper_pos[1] = slack_col[c];
           upper_pos[2] = slack_row[r];
           upper_bound = robust_hinge_min(upper_pos, delta_min);

           int cols_remaining = free_C[j] - c;
           real neutral_logit = -log(cols_remaining);
           this_inv_logit = inv_logit(neutral_logit + lambda[j, r, c]);

           tmp_cell_value[r, c] = lower_bound + this_inv_logit * (upper_bound - lower_bound);
           slack_col[c] = slack_col[c] - tmp_cell_value[r, c];
           slack_row[r] = slack_row[r] - tmp_cell_value[r, c];
           rt = rt - tmp_cell_value[r, c];
           log_det_J += log((upper_bound - lower_bound) * this_inv_logit * (1 - this_inv_logit));
         }
         tmp_cell_value[r, free_C[j]] = slack_row[r];
         rt = rt - tmp_cell_value[r, free_C[j]];
         slack_col[free_C[j]] = slack_col[free_C[j]] - tmp_cell_value[r, free_C[j]];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C[j]];
       }
       for (c in 1:(free_C[j] - 1)){
         tmp_cell_value[free_R[j], c] = slack_col[c];
         rt = rt - tmp_cell_value[free_R[j], c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R[j], c];
         slack_row[free_R[j]] = slack_row[free_R[j]] - tmp_cell_value[free_R[j], c];
       }
       tmp_cell_value[free_R[j], free_C[j]] = rt;

      // =========================================================================
      // 2. DENSE JACOBIAN & LOG-CELL EXTRACTION
      // =========================================================================
      int fr_K = free_R[j] - 1;
      int fc_K = free_C[j] - 1;
      int K_free = fr_K * fc_K;

      // Step A: Calculate the Dense LOR Jacobian on the strictly free cells
      if (K_free > 0) {
        matrix[K_free, K_free] J_lor = rep_matrix(0.0, K_free, K_free);
        real inv_corner = 1.0 / fmax(tmp_cell_value[free_R[j], free_C[j]], 1e-10);

        int row_idx = 1;
        for (r in 1:fr_K) {
          real inv_row_ref = 1.0 / fmax(tmp_cell_value[r, free_C[j]], 1e-10);
          for (c in 1:fc_K) {
            real inv_col_ref = 1.0 / fmax(tmp_cell_value[free_R[j], c], 1e-10);
            real inv_free = 1.0 / fmax(tmp_cell_value[r, c], 1e-10);

            int col_idx = 1;
            for (rp in 1:fr_K) {
              for (cp in 1:fc_K) {
                real derivative = inv_corner;
                if (r == rp && c == cp) derivative += inv_free;
                if (r == rp) derivative += inv_row_ref;
                if (c == cp) derivative += inv_col_ref;
                J_lor[row_idx, col_idx] = derivative;
                col_idx += 1;
              }
            }
            row_idx += 1;
          }
        }
        log_det_J += log_determinant(J_lor);
      }

      // Step B: Build the Log-Scale Matrix to Return
      int fr = 0;
      for(r in 1:R){
        if(row_margins[j, r] > 0){
          fr += 1;
          int fc = 0;
          for(c in 1:C){
            if(col_margins[j, c] > 0){
              fc += 1;
              log_cv_out[j, r, c] = log(fmax(tmp_cell_value[fr, fc], 1e-10));
            } else if(structural_zeros[j, r, c] == 0){
              log_cv_out[j, r, c] = lambda_zero[zero_cell_map[j, r, c]];
            } else {
              log_cv_out[j, r, c] = -200.0;
            }
          }
        } else {
          for(c in 1:C){
            if(structural_zeros[j, r, c] == 0){
              log_cv_out[j, r, c] = lambda_zero[zero_cell_map[j, r, c]];
            } else {
              log_cv_out[j, r, c] = -200.0;
            }
          }
        }
      }
    } // End area loop

    target += log_det_J;
    return log_cv_out;
}

real[,,] ss_assign_ilr_wzeros_hinge_newparam_lp(
    int n_areas, int R, int C,
    matrix row_margins, matrix col_margins,
    real[,,] lambda,
    array[] real lambda_zero,
    array[,,] int zero_cell_map,
    array[,,] int structural_zeros,
    real delta_floor, real delta_min,
    matrix V_ilr) { // ILR Orthonormal Basis Matrix

     // Declare all variables at the top for strict Stan compatibility
     real ILR_jrc[n_areas, R, C - 1]; // <-- Return dimension is now C - 1
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;

     // Variables used during the ILR transformation phase
     int fr;
     int fc;
     vector[C] log_cell_row;
     row_vector[C - 1] ilr_row;
     real active_log_sum;

     lower_pos[1] = 0.0;
     log_det_J = 0;

     // =========================================================================
     // 1. SEQUENTIAL SAMPLING ENGINE (Unchanged from your original logic)
     // =========================================================================
     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);
       free_R = 0;
       free_C = 0;
       for (r in 1:R){
         if(row_margins[j, r] > 0){
           free_R += 1;
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c] > 0){
           free_C += 1;
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }
       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];
       rt = sum(slack_row);
       matrix[free_R, free_C] tmp_cell_value;

       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C - 1)){
           lower_pos[2] = slack_row[r] - sum(tail(slack_col, free_C - c));
           lower_bound = robust_hinge_floor_zero(lower_pos[2], delta_floor);
           upper_pos[1] = slack_col[c];
           upper_pos[2] = slack_row[r];
           upper_bound = robust_hinge_min(upper_pos, delta_min);
           int cols_remaining = free_C - c;
           real neutral_logit = -log(cols_remaining);
           this_inv_logit = inv_logit(neutral_logit + lambda[j, r, c]);
           tmp_cell_value[r, c] = lower_bound + this_inv_logit * (upper_bound - lower_bound);
           slack_col[c] = slack_col[c] - tmp_cell_value[r, c];
           slack_row[r] = slack_row[r] - tmp_cell_value[r, c];
           rt = rt - tmp_cell_value[r, c];
           log_det_J += log((upper_bound - lower_bound) * this_inv_logit * (1 - this_inv_logit));
         }
         tmp_cell_value[r, free_C] = slack_row[r];
         rt = rt - tmp_cell_value[r, free_C];
         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
       }
       for (c in 1:(free_C - 1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt - tmp_cell_value[free_R, c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R, c];
         slack_row[free_R] = slack_row[free_R] - tmp_cell_value[free_R, c];
       }
       tmp_cell_value[free_R, free_C] = rt;

       // =========================================================================
       // 2. REVISED: LATENT ZERO INJECTION & ILR ROTATION
       // =========================================================================
       fr = 0;
       for(r in 1:R){
         if(row_margins[j, r] > 0){
           fr += 1;
           fc = 0;
           active_log_sum = 0.0;

           // Step A: Build the true log-scale composition vector for active rows
           for(c in 1:C){
             if(col_margins[j, c] > 0){
               fc += 1;
               log_cell_row[c] = log(fmax(tmp_cell_value[fr, fc], 1e-10));
               active_log_sum += log_cell_row[c]; // Track only active elements
             } else if(structural_zeros[j, r, c] == 0){
               // Sampling Zero: Inject estimated latent parameter directly
               log_cell_row[c] = lambda_zero[zero_cell_map[j, r, c]];
             } else {
               // Structural Zero
               log_cell_row[c] = -10.0;
             }
           }

           // Step B: Dynamic Jacobian (isolating strictly free cell elements)
           if(fr < free_R){
             log_det_J += log(fmax(row_margins[j, r], 1e-10)) - active_log_sum;
           }

           // Step C: Multiply by the Orthonormal Matrix to map cleanly into ILR Space
           ilr_row = to_row_vector(log_cell_row) * V_ilr;
           for(c in 1:(C - 1)){
             ILR_jrc[j, r, c] = ilr_row[c];
           }

         } else {
           // Step D: Uniformly handle zero-margin rows via ILR to maintain geometry
           for(c in 1:C){
             if(structural_zeros[j, r, c] == 0){
               log_cell_row[c] = lambda_zero[zero_cell_map[j, r, c]];
             } else {
               log_cell_row[c] = -10.0;
             }
           }

           ilr_row = to_row_vector(log_cell_row) * V_ilr;
           for(c in 1:(C - 1)){
             ILR_jrc[j, r, c] = ilr_row[c];
           }
         }
       }
     }

    target += log_det_J; // Ensure Jacobian updates target density internally
    return ILR_jrc;
}


    // constrains using sequential sampling approach described in Chen et. al 2005
    // rather than using sharp bounds, to assist the sampler logistic hinge functions are used to approximate min and floor zero functions
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure.
    // The function then converts the cell_values to Additive Log Ratio components
    // It also adjust the Jacobian to account for the full-transformation from the lambda parameters to the ALR components (which are part of the full log-linear model of the tables).
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     // matrix[n_areas, R] slack_row;
     // matrix[n_areas, C] slack_col;

real[,,] ss_assign_alr_wzeros_hinge_newparam_lp(
    int n_areas, int R, int C,
    matrix row_margins, matrix col_margins,
    real[,,] lambda,
    // array[] real lambda_zero,
    // array[,,] int zero_cell_map,
    array[,,] int structural_zeros,
    real delta_floor, real delta_min){
     real ALR_jrc[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;
     lower_pos[1]=0.0;
     log_det_J = 0;

     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);
       free_R = 0;
       free_C = 0;
       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }
       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];
       rt=sum(slack_row);
       matrix[free_R, free_C] tmp_cell_value;

       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C - 1)){
           lower_pos[2]=slack_row[r]-sum(tail(slack_col, free_C-c));
           lower_bound=robust_hinge_floor_zero(lower_pos[2], delta_floor);
           upper_pos[1]=slack_col[c];
           upper_pos[2]=slack_row[r];
           upper_bound=robust_hinge_min(upper_pos, delta_min);
           int cols_remaining = free_C - c;
           real neutral_logit = -log(cols_remaining);
           this_inv_logit = inv_logit(neutral_logit + lambda[j,r,c]);
           tmp_cell_value[r,c]= lower_bound + this_inv_logit*(upper_bound-lower_bound);
           slack_col[c]=slack_col[c] - tmp_cell_value[r,c];
           slack_row[r]=slack_row[r] - tmp_cell_value[r,c];
           rt = rt - tmp_cell_value[r, c];
           log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
         }
         tmp_cell_value[r, free_C]=slack_row[r];
         rt = rt - tmp_cell_value[r, free_C];
         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
       }
       for (c in 1:(free_C-1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt - tmp_cell_value[free_R, c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R, c];
         slack_row[free_R] = slack_row[free_R] - tmp_cell_value[free_R, c];
       }
       tmp_cell_value[free_R, free_C]=rt;

       // fill ALR_jrc
       int fr = 0;
       for(r in 1:R){
         vector[C] cell_row;
         if(row_margins[j, r] > 0){
           fr += 1;
           int fc = 0;
           for(c in 1:C){
             if(col_margins[j, c] > 0){
               fc += 1;
               cell_row[c] = fmax(tmp_cell_value[fr, fc], 1e-10);
             } else {
               cell_row[c] = 1e-10;
             }
           }
           // ALR Jacobian - only for non-last free rows
           if(fr < free_R){
             log_det_J += log(fmax(row_margins[j, r], 1e-10)) - sum(log(cell_row));
           }
           // ALR transform
           for(c in 1:(C-1)){
             if(col_margins[j, c] > 0){
               ALR_jrc[j, r, c] = log(cell_row[c]) - log(cell_row[C]);
             } else if(structural_zeros[j, r, c] == 0){
               // ALR_jrc[j, r, c] = lambda_zero[zero_cell_map[j, r, c]];
               ALR_jrc[j, r, c] = 0;
             } else {
               ALR_jrc[j, r, c] = 0;
             }
           }
           ALR_jrc[j, r, C] = 0;
         } else {
           // zero margin row
           for(c in 1:(C-1)){
             if(structural_zeros[j, r, c] == 0){
               // ALR_jrc[j, r, c] = lambda_zero[zero_cell_map[j, r, c]];
               ALR_jrc[j, r, c] = 0;
             } else {
               ALR_jrc[j, r, c] = 0;
             }
           }
           ALR_jrc[j, r, C] = 0;
         }
       }
     }

    // target += log_det_J;
    return ALR_jrc;
}



    // constrains using sequential sampling approach described in Chen et. al 2005
    // rather than using sharp bounds, to assist the sampler logistic hinge functions are used to approximate min and floor zero functions
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure.
    // The function then converts the cell_values to Additive Log Ratio components
    // It also adjust the Jacobian to account for the full-transformation from the lambda parameters to the ALR components (which are part of the full log-linear model of the tables).
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     // matrix[n_areas, R] slack_row;
     // matrix[n_areas, C] slack_col;

real[,,] ss_assign_alr_wzeros_hinge_lp(
    int n_areas, int R, int C,
    matrix row_margins, matrix col_margins,
    real[,,] lambda,
    array[] real lambda_zero,
    array[,,] int zero_cell_map,
    array[,,] int structural_zeros,
    real delta_floor, real delta_min){
     real ALR_jrc[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;
     lower_pos[1]=0.0;
     log_det_J = 0;

     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);
       free_R = 0;
       free_C = 0;
       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }
       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];
       rt=sum(slack_row);
       matrix[free_R, free_C] tmp_cell_value;

       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C - 1)){
           lower_pos[2]=slack_row[r]-sum(tail(slack_col, free_C-c));
           lower_bound=robust_hinge_floor_zero(lower_pos[2], delta_floor);
           upper_pos[1]=slack_col[c];
           upper_pos[2]=slack_row[r];
           upper_bound=robust_hinge_min(upper_pos, delta_min);
           this_inv_logit = inv_logit(lambda[j,r,c]);
           tmp_cell_value[r,c]= lower_bound + this_inv_logit*(upper_bound-lower_bound);
           slack_col[c]=slack_col[c] - tmp_cell_value[r,c];
           slack_row[r]=slack_row[r] - tmp_cell_value[r,c];
           rt = rt - tmp_cell_value[r, c];
           log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
         }
         tmp_cell_value[r, free_C]=slack_row[r];
         rt = rt - tmp_cell_value[r, free_C];
         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
       }
       for (c in 1:(free_C-1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt - tmp_cell_value[free_R, c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R, c];
         slack_row[free_R] = slack_row[free_R] - tmp_cell_value[free_R, c];
       }
       tmp_cell_value[free_R, free_C]=rt;

       // fill ALR_jrc
       int fr = 0;
       for(r in 1:R){
         vector[C] cell_row;
         if(row_margins[j, r] > 0){
           fr += 1;
           int fc = 0;
           for(c in 1:C){
             if(col_margins[j, c] > 0){
               fc += 1;
               cell_row[c] = fmax(tmp_cell_value[fr, fc], 1e-10);
             } else {
               cell_row[c] = 1e-10;
             }
           }
           // ALR Jacobian - only for non-last free rows
           if(fr < free_R){
             log_det_J += log(fmax(row_margins[j, r], 1e-10)) - sum(log(cell_row));
           }
           // ALR transform
           for(c in 1:(C-1)){
             if(col_margins[j, c] > 0){
               ALR_jrc[j, r, c] = log(cell_row[c]) - log(cell_row[C]);
             } else if(structural_zeros[j, r, c] == 0){
               ALR_jrc[j, r, c] = lambda_zero[zero_cell_map[j, r, c]];
             } else {
               ALR_jrc[j, r, c] = 0;
             }
           }
           ALR_jrc[j, r, C] = 0;
         } else {
           // zero margin row
           for(c in 1:(C-1)){
             if(structural_zeros[j, r, c] == 0){
               ALR_jrc[j, r, c] = lambda_zero[zero_cell_map[j, r, c]];
             } else {
               ALR_jrc[j, r, c] = 0;
             }
           }
           ALR_jrc[j, r, C] = 0;
         }
       }
     }

    // target += log_det_J;
    return ALR_jrc;
}

real[,,] ss_assign_alr_wzeros_hinge_original_lp (int n_areas, int R, int C, matrix row_margins, matrix col_margins, real[,,] lambda, real delta_floor, real delta_min){
    // constrains using sequential sampling approach described in Chen et. al 2005
    // rather than using sharp bounds, to assist the sampler logistic hinge functions are used to approximate min and floor zero functions
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure.
    // The function then converts the cell_values to Additive Log Ratio components
    // It also adjust the Jacobian to account for the full-transformation from the lambda parameters to the ALR components (which are part of the full log-linear model of the tables).
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     // matrix[n_areas, R] slack_row;
     // matrix[n_areas, C] slack_col;
     real ALR_jrc[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;
     lower_pos[1]=0.0;
     log_det_J = 0;
     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);
       free_R = 0;
       free_C = 0;
       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }
       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];
       rt=sum(slack_row);
       matrix[free_R, free_C] tmp_cell_value;
       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C - 1)){
           lower_pos[2]=slack_row[r]-sum(tail(slack_col, free_C-c));
           lower_bound=robust_hinge_floor_zero(lower_pos[2], delta_floor);
           upper_pos[1]=slack_col[c];
           upper_pos[2]=slack_row[r];
           upper_bound=robust_hinge_min(upper_pos, delta_min);
           this_inv_logit = inv_logit(lambda[j,r,c]);
           tmp_cell_value[r,c]= lower_bound + this_inv_logit*(upper_bound-lower_bound);
           slack_col[c]=slack_col[c] - tmp_cell_value[r,c];
           slack_row[r]=slack_row[r] - tmp_cell_value[r,c];
           rt = rt - tmp_cell_value[r, c];
           // Part 1: Jacobian from lambda to cell_values
           log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
         }
         tmp_cell_value[r, free_C]=slack_row[r];
         rt = rt - tmp_cell_value[r, free_C];
         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
       }
       for (c in 1:(free_C-1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt - tmp_cell_value[free_R, c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R, c];
         slack_row[free_R] = slack_row[free_R] - tmp_cell_value[free_R, c];
       }
       tmp_cell_value[free_R, free_C]=rt;

        // fill ALR_jrc with full R x C dimensions
        int fr = 0;
        for(r in 1:R){
          //compute ALR and Jacobian
          vector[C] cell_row;
          if(row_margins[j, r] > 0){
            fr += 1;
            int fc = 0;
            for(c in 1:C){
              if(col_margins[j, c] > 0){
                fc +=1;
                cell_row[c] = fmax(tmp_cell_value[fr, fc], 1e-10);
                } else {
                  cell_row[c] = 1e-10;
                }
              }
          } else{
            for (c in 1:C){
              cell_row[c] = 1e-10;
            }
          }
          // ALR Jacobian
          if(row_margins[j, r] > 0 && fr < free_R){
            log_det_J += log(fmax(row_margins[j, r], 3e-10)) - sum(log(cell_row));
            }
          // ALR transform
          for(c in 1:(C - 1)){
            ALR_jrc[j, r, c] = log(cell_row[c]) - log(cell_row[C]);
          }
          ALR_jrc[j, r, C] = 0;  // reference column

          }
        }
    target += log_det_J;
    return ALR_jrc;
}



real[,,] ss_assign_log_row_rates_wzeros_hinge_lp (int n_areas, int R, int C, matrix row_margins, matrix col_margins, real[,,] lambda, real delta_floor, real delta_min){
    // constrains using sequential sampling approach described in Chen et. al 2005
    // rather than using sharp bounds, to assist the sampler logistic hinge functions are used to approximate min and floor zero functions
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure.
    // The function then converts the cell_values to log_row_rates
    // It also adjust the Jacobian to account for the full-transformation from the lambda parameters to the log row rates (which are part of the full log-linear model of the tables).
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     // matrix[n_areas, R] slack_row;
     // matrix[n_areas, C] slack_col;
     real log_row_rates[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;
     lower_pos[1]=0.0;
     log_det_J = 0;
     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);
       free_R = 0;
       free_C = 0;
       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }
       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];
       rt=sum(slack_row);
       matrix[free_R, free_C] tmp_cell_value;
       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C - 1)){
           lower_pos[2]=slack_row[r]-sum(tail(slack_col, free_C-c));
           lower_bound=robust_hinge_floor_zero(lower_pos[2], delta_floor);
           upper_pos[1]=slack_col[c];
           upper_pos[2]=slack_row[r];
           upper_bound=robust_hinge_min(upper_pos, delta_min);
           this_inv_logit = inv_logit(lambda[j,r,c]);
           tmp_cell_value[r,c]= lower_bound + this_inv_logit*(upper_bound-lower_bound);
           slack_col[c]=slack_col[c] - tmp_cell_value[r,c];
           slack_row[r]=slack_row[r] - tmp_cell_value[r,c];
           rt = rt - tmp_cell_value[r, c];
           // Part 1: Jacobian from lambda to cell_values
           log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
           // Part 2: Jacobian from cell_values to log row rates
           log_det_J += -log(fmax(tmp_cell_value[r,c], 1e-10));
         }
         tmp_cell_value[r, free_C]=slack_row[r];
         rt = rt - tmp_cell_value[r, free_C];
         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
       }
       for (c in 1:(free_C-1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt - tmp_cell_value[free_R, c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R, c];
         slack_row[free_R] = slack_row[free_R] - tmp_cell_value[free_R, c];
       }
       tmp_cell_value[free_R, free_C]=rt;
       // compute log row rates for all cells
       int fr = 0;
       for (r in 1:R){
         if(row_margins[j, r]>0){
           fr += 1;
           int fc = 0;
           for (c in 1:C){
             if(col_margins[j, c]>0){
               fc += 1;
               log_row_rates[j, r, c] = fmax(
                   log(fmax(tmp_cell_value[fr, fc], 1e-10)) - log(row_margins[j, r]),
                   log(1.0/C) - 2);
             } else {
               log_row_rates[j, r, c] = fmax(
                   log(1e-10) - log(row_margins[j, r]),
                   log(1.0/C) - 2);
             }
           }
         } else {
           // row margin is zero - row rate undefined
           for (c in 1:C){
             log_row_rates[j, r, c] = -200;
           }
         }
       }
     }
    target += log_det_J;
    return log_row_rates;
  }


      real[,,] ss_assign_cvals_wzeros_lp (int n_areas, int R, int C, matrix row_margins, matrix col_margins, real[,,] lambda){
    // constrains using sequential sampling approach described in Chen et. al 2005
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure. The function returns a flattened version of this information, together with an additional value, which is the contains the log determinant of the transform from the unconstrained lambda parameter to the cell values.
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     // matrix[n_areas, R] slack_row;
     // matrix[n_areas, C] slack_col;
     real cell_value[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;

     // slack_row=row_margins;
     // slack_col=col_margins;

     lower_pos[1]=0.0;
     log_det_J = 0;
     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);

       free_R = 0;
       free_C = 0;

       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           // slack_row = append_col(slack_row, row_margins[j, r]);
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           // slack_col = append_col(slack_col, col_margins[j, c]);
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }


       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];

      rt=sum(slack_row);
      matrix[free_R, free_C] tmp_cell_value;

       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C -1 )){
           lower_pos[2]=slack_row[r]-sum(tail(slack_col, free_C-c));
           lower_bound=max(lower_pos);
           upper_pos[1]=slack_col[c];
           upper_pos[2]=slack_row[r];
           upper_bound=min(upper_pos);
           this_inv_logit = inv_logit(lambda[j,r,c]);
           tmp_cell_value[r,c]= lower_bound + this_inv_logit*(upper_bound-lower_bound);
           slack_col[c]=slack_col[c] - tmp_cell_value[r,c];
           slack_row[r]=slack_row[r] - tmp_cell_value[r,c];
           rt = rt - tmp_cell_value[r, c];
           log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
         }
         tmp_cell_value[r, free_C]=slack_row[r];
         rt = rt - tmp_cell_value[r, free_C];
         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
       }

       for (c in 1:(free_C-1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt- tmp_cell_value[free_R, c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R, c];
         slack_row[free_R] = slack_row[free_R] - tmp_cell_value[free_R, c];
       }
       tmp_cell_value[free_R, free_C]=rt;

      int fr = 0;
      for (r in 1:R){
       if(row_margins[j, r]>0){
         fr += 1;
       }
       int fc = 0;
       for (c in 1:C){
         if(col_margins[j, c]>0 && row_margins[j, r]>0){
           fc += 1;
           cell_value[j, r, c] = tmp_cell_value[fr, fc];
         } else {
           cell_value[j, r, c] = 0;
         }
       }

     }
     }

    // Jacobian
    target += log_det_J;
    return cell_value;

  }

      real[,,] ss_assign_cvals_wzeros_hinge_lp (int n_areas, int R, int C, matrix row_margins, matrix col_margins, real[,,] lambda, real delta_floor, real delta_min){
    // constrains using sequential sampling approach described in Chen et. al 2005
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure. The function returns a flattened version of this information, together with an additional value, which is the contains the log determinant of the transform from the unconstrained lambda parameter to the cell values.
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     // matrix[n_areas, R] slack_row;
     // matrix[n_areas, C] slack_col;
     real cell_value[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     int free_R;
     int free_C;
     real this_inv_logit;
     real log_det_J;

     // slack_row=row_margins;
     // slack_col=col_margins;

     lower_pos[1]=0.0;
     log_det_J = 0;
     for (j in 1:n_areas){
       row_vector[R] slack_row_raw = rep_row_vector(0, R);
       row_vector[C] slack_col_raw = rep_row_vector(0, C);

       free_R = 0;
       free_C = 0;

       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           // slack_row = append_col(slack_row, row_margins[j, r]);
           slack_row_raw[free_R] = row_margins[j, r];
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           // slack_col = append_col(slack_col, col_margins[j, c]);
           slack_col_raw[free_C] = col_margins[j, c];
         }
       }


       row_vector[free_R] slack_row = slack_row_raw[1:free_R];
       row_vector[free_C] slack_col = slack_col_raw[1:free_C];

      rt=sum(slack_row);
      matrix[free_R, free_C] tmp_cell_value;

       for (r in 1:(free_R - 1)){
         for (c in 1:(free_C -1 )){
           lower_pos[2]=slack_row[r]-sum(tail(slack_col, free_C-c));
           lower_bound=robust_hinge_floor_zero(lower_pos[2], delta_floor);
           upper_pos[1]=slack_col[c];
           upper_pos[2]=slack_row[r];
           upper_bound=robust_hinge_min(upper_pos, delta_min);
           this_inv_logit = inv_logit(lambda[j,r,c]);
           tmp_cell_value[r,c]= lower_bound + this_inv_logit*(upper_bound-lower_bound);
           slack_col[c]=slack_col[c] - tmp_cell_value[r,c];
           slack_row[r]=slack_row[r] - tmp_cell_value[r,c];
           rt = rt - tmp_cell_value[r, c];
           log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
         }
         tmp_cell_value[r, free_C]=slack_row[r];
         rt = rt - tmp_cell_value[r, free_C];
         slack_col[free_C] = slack_col[free_C] - tmp_cell_value[r, free_C];
         slack_row[r] = slack_row[r] - tmp_cell_value[r, free_C];
       }

       for (c in 1:(free_C-1)){
         tmp_cell_value[free_R, c] = slack_col[c];
         rt = rt- tmp_cell_value[free_R, c];
         slack_col[c] = slack_col[c] - tmp_cell_value[free_R, c];
         slack_row[free_R] = slack_row[free_R] - tmp_cell_value[free_R, c];
       }
       tmp_cell_value[free_R, free_C]=rt;

      int fr = 0;
      for (r in 1:R){
       if(row_margins[j, r]>0){
         fr += 1;
       }
       int fc = 0;
       for (c in 1:C){
         if(col_margins[j, c]>0 && row_margins[j, r]>0){
           fc += 1;
           cell_value[j, r, c] = tmp_cell_value[fr, fc];
         } else {
           cell_value[j, r, c] = 0;
         }
       }

     }
     }

    // Jacobian
    target += log_det_J;
    return cell_value;

  }

      real[,,] ss_log_cvals_wzeros_lp (int n_areas, int R, int C, matrix row_margins, matrix col_margins, real[,,] cell_values){

     real log_cell_values[n_areas, R, C];
     int free_R[n_areas];
     int free_C[n_areas];
     real log_det_J = 0;

     for (j in 1:n_areas){
       free_R[j] = 0;
       free_C[j] = 0;

       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R[j] += 1;
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C[j] += 1;
         }
       }


      int fr = 0;
      for (r in 1:R){
       if(row_margins[j, r]>0){
         fr += 1;
       }
       int fc = 0;
       for (c in 1:C){
         if(col_margins[j, c]>0 && row_margins[j, r]>0){
           fc += 1;
           log_cell_values[j, r, c] = log(cell_values[j, r, c]);
           if(fc<free_R[c] && fr<free_R[j]){
             log_det_J += -log_cell_values[j, r, c];
           }
         } else {
           log_cell_values[j, r, c] = log(.001);
         }
       }

     }
     }

    // Jacobian
    target += log_det_J;
    return log_cell_values;

  }
