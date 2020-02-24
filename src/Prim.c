// -*-  mode: c    ; coding: koi8   -*- ----------------------------------------

//------------------------------------------------------------------------------

#include <stdlib.h>   
#include <string.h>

#include "gsl_multimin.h"
#include "gsl_blas.h"  // обязательно, а то идут "nun"                               
#include "gsl_vector.h"

#include "XXX.h"

//------------------------------------------------------------------------------

//The following example function defines a simple paraboloid with two
//parameters,

//------------------------------------------------------------------------------
/*  Paraboloid centered on (dp[0],dp[1])                                      */
//------------------------------------------------------------------------------
double
my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double *dp = (double *) params; // координаты цeнтра параболы
  
  x = gsl_vector_get (v, 0);
  y = gsl_vector_get (v, 1);
 
  //fprintf (stderr, "my_f : x= %f  y= %f \n", x, y);
  //fprintf (stderr, "my_f : dp_0= %f  dp_1= %f \n", dp[0], dp[1]);

  return (10.0 * (x - dp[0]) * (x - dp[0]) +
          20.0 * (y - dp[1]) * (y - dp[1]) + 30.0); 
}
//------------------------------------------------------------------------------
/*  The gradient of f, df = (df/dx, df/dy).                                   */
//------------------------------------------------------------------------------
void 
my_df (const gsl_vector *v, void *params, 
       gsl_vector *df)
{
  double x, y;
  double *dp = (double *) params; // координаты цeнтра параболы ????? !!!!
  
  x = gsl_vector_get (v, 0);
  y = gsl_vector_get (v, 1);
 
  //fprintf (stderr, "my_df: x= %f  y= %f \n", x, y);
  //fprintf (stderr, "my_df: dp_0= %f  dp_1= %f \n", dp[0], dp[1]);

  // ищeм градиeнты:

  gsl_vector_set (df, 0, 20.0 * (x - dp[0]));
  gsl_vector_set (df, 1, 40.0 * (y - dp[1]));
}
//------------------------------------------------------------------------------
/* Compute both f and df together. */
//------------------------------------------------------------------------------
void 
my_fdf (const gsl_vector *x, void *params, 
        double *f, gsl_vector *df) 
{

  *f = my_f (x, params); 

  my_df (x, params, df);

/*   fprintf (stderr, "my_fdf:  f  = %f \n", *f ); */
/*   fprintf (stderr, "my_fdf: df_0  = %f \n", gsl_vector_get (df, 0) ); */
/*   fprintf (stderr, "my_fdf: df_1  = %f \n", gsl_vector_get (df, 1) ); */
}
//------------------------------------------------------------------------------
int
main (int argc, char **argv)
{
  int iter_num = 100;
  //int iter_num = 1;

  char *name_minimizer;
  int   iter = 0;
  int   status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_multimin_function_fdf  my_func;
  my_func.f   = &my_f;
  my_func.df  = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n   = 2;

  double par[2] = { 1.0, 2.0 };  // координаты цeнтра параболы, он жe минимум
  my_func.params = &par;

  /* Starting point, x = (5, 7) */
  gsl_vector *x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 5.0);
  gsl_vector_set (x, 1, 7.0);

  if (argc == 1) {
    fprintf (stderr, "\n");
    fprintf (stderr, "Prim step|pr|fr|bfgs \n");
    fprintf (stderr, "\n");
    exit (1);
  }

  name_minimizer = argv[1];
  if      (! strcmp (name_minimizer, "step"))  T = gsl_multimin_fdfminimizer_steepest_descent;
  else if (! strcmp (name_minimizer, "pr"))    T = gsl_multimin_fdfminimizer_conjugate_pr;
  else if (! strcmp (name_minimizer, "fr"))    T = gsl_multimin_fdfminimizer_conjugate_fr;
  else if (! strcmp (name_minimizer, "bfgs"))  T = gsl_multimin_fdfminimizer_vector_bfgs;

  s = gsl_multimin_fdfminimizer_alloc (T, 2);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

  //-------------------------------------------
  do {
    iter++;

    status = gsl_multimin_fdfminimizer_iterate (s);    
    if (status)
      break;
    
    status = gsl_multimin_test_gradient (s->gradient, 1e-3);
    
    if (status == GSL_SUCCESS)
      printf ("Minimum found at:\n");
    
    printf ("%5d)  ", iter);
    xxx_vector_printf_line ("x= ", s->x);
    printf ("%10.5f \n", s->f);    

  }
  while (status == GSL_CONTINUE && iter < iter_num /* 3  *//* 100 */);
  //-------------------------------------------
  
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return 0;
}
//------------------------------------------------------------------------------
//The initial step-size is chosen as 0.01, a conservative estimate in this
//case, and the line minimization parameter is set at 0.0001.  The program
//terminates when the norm of the gradient has been reduced below
//0.001. The output of the program is shown below,

//         x       y         f
//    1 4.99629 6.99072  687.84780
//    2 4.98886 6.97215  683.55456
//    3 4.97400 6.93501  675.01278
//    4 4.94429 6.86073  658.10798
//    5 4.88487 6.71217  625.01340
//    6 4.76602 6.41506  561.68440
//    7 4.52833 5.82083  446.46694
//    8 4.05295 4.63238  261.79422
//    9 3.10219 2.25548   75.49762
//   10 2.85185 1.62963   67.03704
//   11 2.19088 1.76182   45.31640
//   12 0.86892 2.02622   30.18555
//Minimum found at:
//   13 1.00000 2.00000   30.00000

//------------------------------------------------------------------------------
