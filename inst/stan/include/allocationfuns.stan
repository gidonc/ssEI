  real[,,] ss_assign_cvals_lp (int n_areas, int R, int C, matrix row_margins, matrix col_margins, real[,,] lambda){
    // constrains using sequential sampling approach described in Chen et. al 2005
    // function to transform unconstrained (R-1)*(C-1) unconstrained parameters (lambda) into an RxC matrix with fixed row and column margins. The function completes the internal structure  across matrices from n_areas regions
    // The function completes cell_value (which is n_area R*C matricies) because indexing errors are easier to spot with this structure. The function returns a flattened version of this information, together with an additional value, which is the contains the log determinant of the transform from the unconstrained lambda parameter to the cell values.
    // n_areas is the number of matrices which need the internal structure completing
    // R is the number of rows
    // C is the number of columns
    // row_margins is a matrix of the row margins in each area
    // col_margins is the matrix of the column margins in each area
     matrix[n_areas, R] slack_row;
     matrix[n_areas, C] slack_col;
     real cell_value[n_areas, R, C];
     vector[2] lower_pos;
     vector[2] upper_pos;
     real lower_bound;
     real upper_bound;
     real rt;
     real this_inv_logit;
     real log_det_J;

     slack_row=row_margins;
     slack_col=col_margins;

     lower_pos[1]=0.0;
     log_det_J = 0;
     for (j in 1:n_areas){
       rt=sum(row(slack_row, j));

       for (r in 1:(R-1)){
         for (c in 1:(C-1)){
           lower_pos[2]=slack_row[j, r]-sum(tail(row(slack_col, j), C-c));
           lower_bound=max(lower_pos);
           upper_pos[1]=slack_col[j, c];
           upper_pos[2]=slack_row[j, r];
           upper_bound=min(upper_pos);
           this_inv_logit = inv_logit(lambda[j,r,c]);
           cell_value[j, r,c]= lower_bound + this_inv_logit*(upper_bound-lower_bound);
           slack_col[j, c]=slack_col[j, c]-cell_value[j, r,c];
           slack_row[j, r]=slack_row[j, r]-cell_value[j, r,c];
           rt=rt-cell_value[j, r,c];
           log_det_J += log((upper_bound - lower_bound)*this_inv_logit*(1-this_inv_logit));
         }
         cell_value[j, r, C]=slack_row[j, r];
         rt=rt-cell_value[j, r, C];
         slack_col[j, C]=slack_col[j, C]-cell_value[j, r, C];
         slack_row[j, r]=slack_row[j, r]-cell_value[j, r, C];
       }

       for (c in 1:(C-1)){
         cell_value[j, R, c]=slack_col[j, c];
         rt=rt-cell_value[j, R, c];
         slack_col[j, c]=slack_col[j, c]-cell_value[j, R, c];
         slack_row[j, R]=slack_row[j, R]-cell_value[j, R, c];
       }
       cell_value[j, R, C]=rt;
     }
     print(cell_value);

    // Jacobian
    target += log_det_J;
    return cell_value;
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
       row_vector[0] slack_row;
       row_vector[0] slack_col;
       free_R = 0;
       free_C = 0;
       for (r in 1:R){
         if(row_margins[j, r]>0){
           free_R += 1;
           slack_row = append_col(slack_row, row_margins[j, r]);
         }
       }
       for (c in 1:C){
         if(col_margins[j, c]>0){
           free_C += 1;
           slack_col = append_col(slack_col, col_margins[j, c]);
         }

       }

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

