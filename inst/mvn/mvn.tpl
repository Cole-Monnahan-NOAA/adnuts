GLOBALS_SECTION
  #include <mvnorm.h>
  MVNORM_t neg_log_density;

DATA_SECTION
  init_number dim;
  init_matrix covar(1,dim,1,dim);

PARAMETER_SECTION
  init_vector mu(1,dim);

objective_function_value jnll;

PRELIMINARY_CALCS_SECTION
  neg_log_density.setSigma(covar);

PROCEDURE_SECTION
  dvector x(1,dim);
  for(int i=1; i<=dim; i++) x(i)=0;
  jnll=neg_log_density(x-mu);

TOP_OF_MAIN_SECTION
  arrmblsize=20000000;
