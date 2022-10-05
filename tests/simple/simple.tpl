// modified "simple" example packaged with ADMB, slight modificationst to work as test model for MCMC

DATA_SECTION
 !!CLASS ofstream report1("mceval.dat")

  init_int nobs
  init_vector Y(1,nobs)
  init_vector x(1,nobs)
PARAMETER_SECTION
  init_bounded_number a(-10,8);
  init_bounded_number b(-19,15);
  vector pred_Y(1,nobs)
  sdreport_number aa
  objective_function_value f
PROCEDURE_SECTION
 aa=a;
  pred_Y=a*x+b;
  f=(norm2(pred_Y-Y));
  f=nobs/2.*log(f);    // make it a likelihood function so that
                       // covariance matrix is correct
 if (mceval_phase()) report1 << a << "," << b << "," << f <<endl;
