//------------------------------------------------------------------------------
/* multimin/vector_bfgs.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 */

/* vector_bfgs.c -- Limited memory Broyden-Fletcher-Goldfarb-Shanno
   conjugate gradient method */

/* Modified by Brian Gough to use single iteration structure */
//------------------------------------------------------------------------------

#include "config.h"

#include "gsl_multimin.h"
#include "gsl_blas.h"

#include "XXX.h"
#include "m_directional_minimize.h" 


//------------------------------------------------------------------------------
void
direction_update_bfgs (xxx_conjugate_state_t *state,
                       gsl_vector *x, 
                       gsl_vector *p, 
                       gsl_vector *gradient,
                       double g1norm, double g0norm, gsl_vector *g0,
                       gsl_vector  *x0
                       )
{

  // This is the BFGS update: 
  // p' = g1 - A dx - B dg 
  // A  = - (1+ dg.dg/dx.dg) B + dg.g/dx.dg 
  // B  = dx.g/dx.dg 

  gsl_vector *dx0 = state->dx0;
  gsl_vector *dg0 = state->dg0;

  double dxg, dgg, dxdg, dgnorm, A, B;

  //fprintf (stderr, "direction_update_bfgs...  \n");

  // dx0 = x - x0 
  gsl_vector_memcpy (dx0, x);
  gsl_blas_daxpy (-1.0, x0, dx0);

/*   xxx_vector_fprintf (stderr, "x  =   ", x);  */
/*   xxx_vector_fprintf (stderr, "x0 =   ", x0);  */
/*   xxx_vector_fprintf (stderr, "dx0=   ", dx0);  */

  // dg0 = g - g0 
  gsl_vector_memcpy (dg0, gradient);
  gsl_blas_daxpy (-1.0, g0, dg0);

  gsl_blas_ddot (dx0, gradient, &dxg);
  gsl_blas_ddot (dg0, gradient, &dgg);
  gsl_blas_ddot (dx0, dg0, &dxdg);

  dgnorm = gsl_blas_dnrm2 (dg0);

/*   xxx_vector_fprintf (stderr, "dg0=   ", dg0);  */
/*              fprintf (stderr, "dxdg=  %g \n", dxdg); */

  if (dxdg != 0) 
  {
    //fprintf (stderr, "1.. \n") ;
    B = dxg / dxdg;
    A = -(1.0 + dgnorm * dgnorm / dxdg) * B + dgg / dxdg;
  }
  else
  {
    //fprintf (stderr, "2.. \n") ;
    B = 0;
    A = 0; 
  }

  gsl_vector_memcpy (p, gradient);
  gsl_blas_daxpy (-A, dx0, p);
  gsl_blas_daxpy (-B, dg0, p);

  state->pnorm = gsl_blas_dnrm2 (p);

/*   xxx_vector_fprintf (stderr, "p=     ", p);  */
/*              fprintf (stderr, "pnorm= %g \n", state->pnorm); */

  return;
}
//------------------------------------------------------------------------------
void
choose_new_direction_bfgs (xxx_conjugate_state_t *state,
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
                        direction_update_bfgs // !!
                        );

  //
  gsl_vector_memcpy (x0, x);
  state->g0norm = gsl_blas_dnrm2 (g0);

}
//------------------------------------------------------------------------------
int
vector_bfgs_alloc (void *vstate, size_t n)
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  //---------------------------------------------------------
  xxx_conjugate_alloc (vstate, n);
  //---------------------------------------------------------

  state->x0 = gsl_vector_calloc (n);
  if (state->x0 == 0) {
    GSL_ERROR ("failed to allocate space for g0", GSL_ENOMEM);
  }

  state->dx0 = gsl_vector_calloc (n);
  if (state->dx0 == 0) {
    GSL_ERROR ("failed to allocate space for g0", GSL_ENOMEM);
  }
  
  state->dg0 = gsl_vector_calloc (n);
  if (state->dg0 == 0) {
    GSL_ERROR ("failed to allocate space for g0", GSL_ENOMEM);
  }

  return GSL_SUCCESS;
}
//------------------------------------------------------------------------------
int
vector_bfgs_set  (void *vstate, 
                 gsl_multimin_function_fdf *fdf,
                 const gsl_vector *x, 
                 double     *f, 
                 gsl_vector *gradient,
                 double step_size, double tol)
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  xxx_conjugate_set (vstate, 
                     fdf,
                     x, 
                     f, 
                     gradient,
                     step_size, tol);


  gsl_vector_memcpy (state->x0, x);

/*   fprintf (stderr, "vector_bfgs_set...  \n"); */
/*   xxx_vector_fprintf (stderr, "x0 =   ", state->x0);  */

  return GSL_SUCCESS;
}
//------------------------------------------------------------------------------
void
vector_bfgs_free (void *vstate)
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  xxx_conjugate_free (vstate);

  gsl_vector_free (state->dg0);
  gsl_vector_free (state->dx0);
  gsl_vector_free (state->x0);
}
//------------------------------------------------------------------------------
int
vector_bfgs_iterate (void *vstate, 
                     gsl_multimin_function_fdf * fdf,
                     gsl_vector *x, double *f,
                     gsl_vector *gradient, gsl_vector *dx)
{

  return (
  xxx_conjugate_iterate (vstate, 
                         fdf,
                         x, f,
                         gradient, dx,
                         choose_new_direction_bfgs
                         ));
}
//------------------------------------------------------------------------------


static const gsl_multimin_fdfminimizer_type vector_bfgs_type = {
  "vector_bfgs",                /* name */
  sizeof (xxx_conjugate_state_t),

  &vector_bfgs_alloc,
  &vector_bfgs_set,
  &vector_bfgs_iterate,
  &xxx_conjugate_restart,
  &vector_bfgs_free
};

const gsl_multimin_fdfminimizer_type
  * gsl_multimin_fdfminimizer_vector_bfgs = &vector_bfgs_type;

//------------------------------------------------------------------------------
