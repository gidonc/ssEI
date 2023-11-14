
  matrix[R,C] row_rate[n_areas];
  matrix[R, C] overall_cell_values;
  matrix[R, C] overall_row_rate;

  for (j in 1:n_areas){
    for(r in 1:R){
      if(row_margins[j, r] == 0){
        // avoid division by zero erors
        row_rate[j, r] = rep_row_vector(-1, C);
      } else {
        row_rate[j, r] = to_row_vector(cell_values[j, r]) / row_margins[j, r];
      }
    }
  }
  for(r in 1:R){
    for(c in 1:C){
      overall_cell_values[r, c] = sum(cell_values[1:n_areas, r, c]);
    }
  }
  for(r in 1:R){
    real this_row_sum = sum(overall_cell_values[r]);
    if(this_row_sum == 0){
      overall_row_rate[r] = rep_row_vector(-1, C);
    } else {
      overall_row_rate[r] = overall_cell_values[r] / this_row_sum;
    }
  }
