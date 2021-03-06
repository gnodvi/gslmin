// -*-  mode: c    ; coding: koi8   -*- ----------------------------------------

//------------------------------------------------------------------------------

#include <stdlib.h>
#include <string.h>

#include "gsl_multimin.h"

//------------------------------------------------------------------------------

static gsl_multimin_function_fdf  my_func;

static double (* global_my_f)  (const gsl_vector * x, void * params);
static   void (* global_my_df) (const gsl_vector * x, void * params, gsl_vector * df);

//------------------------------------------------------------------------------
/* Compute both f and df together. */
//------------------------------------------------------------------------------
void 
my_fdf (const gsl_vector *x, void *params, 
        double *f, gsl_vector *df) 
{

  //*f = my_f (x, params); 
  // my_df (x, params, df);

  *f = global_my_f (x, params); 
  global_my_df (x, params, df);

/*   fprintf (stderr, "my_fdf:  f  = %f \n", *f ); */
/*   fprintf (stderr, "my_fdf: df_0  = %f \n", gsl_vector_get (df, 0) ); */
/*   fprintf (stderr, "my_fdf: df_1  = %f \n", gsl_vector_get (df, 1) ); */
  // �����e ��� � ������� �e��e �������e��� ������ ���� ���, ������� �������e��e
}
//------------------------------------------------------------------------------
// ����� �e ��������� ��� ��������� "gsl_multimin_fdfminimizer_type" �e�����
//------------------------------------------------------------------------------
gsl_multimin_fdfminimizer *
my_gsl_multimin_fdfminimizer_alloc (char *type, size_t n)
{
  gsl_multimin_fdfminimizer *s;
  const gsl_multimin_fdfminimizer_type *T;

  if      (!strcmp (type, "gsl_multimin_fdfminimizer_conjugate_fr"))
    T = gsl_multimin_fdfminimizer_conjugate_fr;

  else if (!strcmp (type, "gsl_multimin_fdfminimizer_conjugate_pr"))
    T = gsl_multimin_fdfminimizer_conjugate_pr;

  else if (!strcmp (type, "gsl_multimin_fdfminimizer_steepest_descent"))
    T = gsl_multimin_fdfminimizer_steepest_descent;

  else if (!strcmp (type, "gsl_multimin_fdfminimizer_vector_bfgs"))
    T = gsl_multimin_fdfminimizer_vector_bfgs;

  s = gsl_multimin_fdfminimizer_alloc (T, /* 2 */n);

  return (s);
}
//------------------------------------------------------------------------------
void
my_gsl_multimin_fdfminimizer_set (
                   gsl_multimin_fdfminimizer * s,

                   double (* my_f)  (const gsl_vector * x, void * params),
                   void   (* my_df) (const gsl_vector * x, void * params, gsl_vector * df),

                                   const gsl_vector * x,
                                   const gsl_vector * par,
                                   double step_size, double tol
                                  )
{

  my_func.f   = global_my_f  = my_f;
  my_func.df  = global_my_df = my_df; 

  my_func.fdf = &my_fdf;

  my_func.n   = /* 2 */x->size;

  my_func.params = (void *) par;

  //return (
  gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tol);
  //         );
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
