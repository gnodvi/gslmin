//------------------------------------------------------------------------------
/* multimin/conjugate_pr.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 */


/* conjugate_pr.c -- Conjugate gradient Polak-Ribiere algorithm */

/* Modified by Brian Gough to use single iteration structure */
//------------------------------------------------------------------------------

#include "config.h"

#include "gsl_multimin.h"
#include "gsl_blas.h"

#include "XXX.h"
#include "m_directional_minimize.h" 
#include "m_conjugate_pr.h"

//------------------------------------------------------------------------------
double
calc_beta_pr (double g0norm, double g1norm,
              gsl_vector *gradient, gsl_vector *g0)
{

  double g0g1, beta;

  gsl_blas_daxpy (-1.0, gradient, g0); // g0'  = g0 - g1 
  gsl_blas_ddot(g0, gradient, &g0g1);  // g1g0 = (g0-g1).g1 
  beta = g0g1 / (g0norm*g0norm);       // beta = -((g1 - g0).g1)/(g0.g0) 

  return (beta);
}
//------------------------------------------------------------------------------
void
direction_update_pr (xxx_conjugate_state_t *state,

                     gsl_vector *x, 
                     gsl_vector *p, 
                     gsl_vector *gradient,
                     double g1norm, double g0norm, gsl_vector *g0,
                     gsl_vector  *x0
                     )
{

  //  p' = g1 - beta * p 
  double beta = calc_beta_pr (g0norm, g1norm, gradient, g0);
  calc_p (state, beta, p, gradient);

}
//------------------------------------------------------------------------------
void
choose_new_direction_pr (xxx_conjugate_state_t *state,
                         gsl_vector *x, 
                         gsl_vector *p, 
                         gsl_vector *gradient,
                         double g1norm, double g0norm, gsl_vector *g0,
                         gsl_vector  *x0
                         )
{
  
  choose_new_direction (state,
                        x, 
                        p, 
                        gradient,
                        g1norm, g0norm, g0,
                        x0,
                        direction_update_pr // !!
                        );
  //
  state->g0norm = g1norm;

}
//------------------------------------------------------------------------------
int
conjugate_pr_iterate (void *vstate, 
                      gsl_multimin_function_fdf *fdf,
                      gsl_vector *x, double *f,
                      gsl_vector *gradient, gsl_vector *dx)
{

  return (
  xxx_conjugate_iterate (vstate, 
                         fdf,
                         x, f,
                         gradient, dx,
                         choose_new_direction_pr
                         ));
}
//------------------------------------------------------------------------------


static const gsl_multimin_fdfminimizer_type conjugate_pr_type = {
  "conjugate_pr",               /* name */
  sizeof (xxx_conjugate_state_t),
  &xxx_conjugate_alloc,
  &xxx_conjugate_set,
  &conjugate_pr_iterate,
  &xxx_conjugate_restart,
  &xxx_conjugate_free
};

const gsl_multimin_fdfminimizer_type
  * gsl_multimin_fdfminimizer_conjugate_pr = &conjugate_pr_type;

//------------------------------------------------------------------------------
