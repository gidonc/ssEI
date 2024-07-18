
matrix lpivot (matrix m, int pivrow, int pivcol){
  matrix[rows(m), cols(m)] mout;
  real pivval = m[pivrow, pivcol];
  real pivremval;

  mout = m;

  mout[pivrow, 1:cols(m)] = mout[pivrow, 1:cols(m)]./pivval;

  for (r in 1:rows(m)){
    if (r != pivrow){
      pivremval = mout[r, pivcol];
      mout[r] = mout[r] - (pivremval .* mout[pivrow]);
    }
  }
  return(mout);
}

matrix addart(matrix m){

  int norig = cols(m);
  int rm = rows(m);
  matrix[rm + 1, norig + rm] mout;
  matrix[rm, rm] art = diag_matrix(rep_vector(1, rm));
  row_vector[norig + rm] objfun;

  objfun = append_col(append_col(rep_row_vector(0, norig - 1), rep_row_vector(1, rm)), 0);

  // add row with phase 1 objective function
  mout = append_row(append_col(append_col(m[1:rm, 1:(norig - 1)], art),  m[1:rm, norig]), objfun);
  for (r in 1:(rows(mout) - 1)){
    mout[rows(mout)] = mout[rows(mout)] - mout[r];
  }

  return(mout);
}

matrix dropart(matrix m){

  int rm = rows(m) - 1;
  int cm = cols(m) - rm;
  matrix[rm, cm] mout;


  // final row is objective function so drop that
  // the n non-artifical variables are in colums 1:(n - 1) and the final column

  mout = append_col(m[1:rm, 1:(cm - 1)], m[1:rm, cols(m)]);
  return(mout);
}

int mrtidx(vector avec, vector bvec){
  // returns the index of the row which meets the minimum ratio test
  // if there are only negative values it returns 0 (this would mean that the problem is unbounded)
  // if there are multiple equal max negative values it returns the index of the first occurance n the row_vector
  int idx = 0;
  real curmrt = -1;
  real thismrt;
  for (i in 1:num_elements(avec)){
    if(avec[i] != 0){
      thismrt = bvec[i]/avec[i];
      if((thismrt < curmrt || curmrt == -1)  && thismrt>0){
        curmrt = thismrt;
        idx = i;
      }
    }
  }
  return(idx);
}

int minnegidx(row_vector row_vec){
  // returns the index of the minimum negative value of a row_vector (excluding zero)
  // if there are no negative values it returns 0
  // if there are multiple equal minimum negative values it returns the index of the first occurance n the row_vector
  int idx = 0;
  real curmin = 0;
  for (i in 1:num_elements(row_vec)){
    if(row_vec[i] < curmin){
      curmin = row_vec[i];
      idx = i;
    }
  }
  return(idx);
}

matrix run_simplex_alg(matrix m){
  int rpiv;
  int cpiv;
  matrix[rows(m), cols(m)] mout;

  cpiv = minnegidx(row(m, rows(m)));
  mout = m;

  while(minnegidx(mout[rows(mout), 1:(cols(mout) - 1)]) > 0){
    cpiv = minnegidx(mout[rows(mout), 1:(cols(mout) - 1)]);
    rpiv = mrtidx(mout[1:(rows(mout) - 1), cpiv], mout[1:(rows(mout) - 1), cols(mout)]);
    if(rpiv == 0) {
      reject("minimum ratio test found no positive values. Indicates an illegal unbounded maximization problem. Current matrix values are :", m);}
    mout = lpivot(mout, rpiv, cpiv);
  }
  return(mout);
}

int which_basic(vector vec){
  // function to find the index of the row that gives the value of a column which is in the basis of a linear programming solution
  // returns 0 if the column does not belong to a variable which is in the basis
  int not_zero[num_elements(vec)];
  int is_basic;
  int idx;

  for(i in 1:num_elements(vec)){
    if(vec[i] != 0){
      not_zero[i] = 1;
    } else{
      not_zero[i] = 0;
    }
  }

  idx = 0;
  is_basic = sum(not_zero) == 1;
  if(is_basic == 1){
    for(i in 1:num_elements(not_zero)){
      if(not_zero[i] != 0){
        idx = i;
      }
    }
  }
  return(idx);
}

vector get_var_lims(matrix m, int cm){
  // find limits for one variable from a matrix m giving a feasible starting point. Uses simplex algorithm with objective fuction which contains just one variable where the function to be maximized in f(x) =x (gives max) and f(x) = -x (gives min)
  int nvars = cols(m) - 1;
  int is_basic;
  int not_zero;
  int use_row;
  vector[2] var_lims;
  matrix[1,1] mat_min;
  matrix[1,1] mat_max;
  row_vector[nvars + 1] obj_max = rep_row_vector(0, nvars + 1);
  row_vector[nvars + 1] obj_min = rep_row_vector(0, nvars + 1);

  obj_max[cm] = -1;
  obj_min[cm] = 1;

  use_row = which_basic(col(m, cm));
  if(use_row > 0){
    obj_max = obj_max + m[use_row];
    obj_min = obj_min - m[use_row];
  }
  mat_min = run_simplex_alg(append_row(m, obj_min));
  mat_min = run_simplex_alg(append_row(m, obj_max));

  use_row = which_basic(col(mat_min, cm));
  if(use_row >0){
    var_lims[1] = mat_min[use_row, cols(mat_min)]/mat_min[use_row, cm];
  } else{
    var_lims[1] = 0;
  }

  use_row = which_basic(col(mat_max, cm));
  if(use_row >0){
    var_lims[2] = mat_max[use_row, cols(mat_max)]/mat_max[use_row, cm];
  } else{
    var_lims[1] = 0;
  }

  return(var_lims);
}
