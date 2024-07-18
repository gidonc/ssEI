
functions{
  #include include/allocationfuns.stan
  #include include/realpdf.stan
  #include include/lkjonionfun.stan
  #include include/linprogfun.stan

}
data{
 int<lower=0> J;
 int<lower=0> R;  // number of rows
 int<lower=0> C;  // number of columns
 int<lower=0, upper=1> lflag_has_target;
 matrix<lower=0>[J, R] row_margins; // the row margins in each area
 matrix<lower=0>[J, C] col_margins; // the column margins in each area
 matrix<lower=0>[R, C] overall_tots; // the overall total target distribution (if known - used only if lflag_has_target==1)

}
transformed data{
  int nconst = J*R + J*C + R*C;
  int ncells = J*R*C;
  array[J, R, C] int initvarpos;
  matrix[nconst, ncells + 1] cequns = rep_matrix(0, nconst, ncells + 1);
  matrix[nconst + 1, ncells + nconst + 1] cequnsart = rep_matrix(0, nconst + 1, ncells + nconst + 1);
  int curvarpos = 0;
  int the_row = 0;
  for (j in 1:J){
    for (r in 1:R){
      for(c in 1:C){
        curvarpos += 1;
        initvarpos[j, r, c] = curvarpos;
      }
    }
  }

  for(j in 1:J){
    for(r in 1:R){
      the_row += 1;
      cequns[the_row, ncells + 1] = row_margins[j, r];
      for(c in 1:C){
        cequns[the_row, initvarpos[j, r, c]] = 1;
      }
    }
  }
  for(j in 1:J){
    for(c in 1:C){
      the_row += 1;
      cequns[the_row, ncells + 1] = col_margins[j, c];
      for(r in 1:R){
        cequns[the_row, initvarpos[j, r, c]] = 1;
      }
    }
  }
  for(r in 1:R){
    for(c in 1:C){
      the_row += 1;
      cequns[the_row, ncells + 1] = overall_tots[r, c];
      for(j in 1:J){
        cequns[the_row, initvarpos[j, r, c]] = 1;
      }
    }
  }

  cequnsart = addart(cequns);
  cequnsart = run_simplex_alg(cequnsart);
  cequns = dropart(cequnsart);



}
parameters{
  real w;
}
transformed parameters{

}

model{
  w ~ normal(0, 1);
}
generated quantities{


}
