// -*-  mode: c    ; coding: koi8   -*- ----------------------------------------

//------------------------------------------------------------------------------

#include "config.h"

#include "gsl_multimin.h"
#include "gsl_blas.h"

#include "XXX.h"
#include "m_directional_minimize.h" 

int _DEBUG = 0;
int D_PRINT = 0/* 1 */;

//------------------------------------------------------------------------------
void
xxx_vector_fprintf_line (FILE *stream, char *name, gsl_vector *v)
{
  int i, n = v->size ;

  fprintf (stream, "%s(", name);

  for (i = 0; i < n; i++) {
    fprintf (stream, "% f", (v->data)[i]);
    if (i != n-1)
    fprintf (stream, " ");
  }

  fprintf (stream, ")");

}
//------------------------------------------------------------------------------
void
xxx_vector_printf_line (char *name, gsl_vector *v)
{

  xxx_vector_fprintf_line (stdout, name, v);

}
//------------------------------------------------------------------------------
void
xxx_vector_fprintf (FILE *stream, char *name, gsl_vector *v)
{
  xxx_vector_fprintf_line (stream, name, v);
  fprintf (stream, " \n");

}
//------------------------------------------------------------------------------
void
xxx_vector_printf (char *name, gsl_vector *v)
{

  xxx_vector_fprintf (stdout, name, v);

}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void
take_step (const gsl_vector *x, const gsl_vector *p,
           double step, double lambda, 
           gsl_vector *x1, gsl_vector *dx)
{
  double alpha = -step * lambda;
  //fprintf (stderr, "take_step:  step= %f  lambda= %f \n", step, lambda);

  gsl_vector_set_zero (dx);      // dx = 0 
  gsl_blas_daxpy (alpha, p, dx); // dx = dx + alpha * p

  gsl_vector_memcpy (x1, x);     // x1 = x
  gsl_blas_daxpy (1.0, dx, x1);  // x1 = x1 + 1.0 * dx

}
//------------------------------------------------------------------------------
int
xxx_conjugate_alloc (void *vstate, size_t n)
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  state->x1 = gsl_vector_calloc (n);
  if (state->x1 == 0) {
    GSL_ERROR ("failed to allocate space for x1", GSL_ENOMEM);
  }
  
  state->dx1 = gsl_vector_calloc (n);
  if (state->dx1 == 0) {
    GSL_ERROR ("failed to allocate space for dx1", GSL_ENOMEM);
  }
  
  state->x2 = gsl_vector_calloc (n);
  if (state->x2 == 0) {
    GSL_ERROR ("failed to allocate space for x2", GSL_ENOMEM);
  }
  
  state->p = gsl_vector_calloc (n);
  if (state->p == 0) {
    GSL_ERROR ("failed to allocate space for p", GSL_ENOMEM);
  }
  
  state->g0 = gsl_vector_calloc (n);
  if (state->g0 == 0) {
    GSL_ERROR ("failed to allocate space for g0", GSL_ENOMEM);
  }

  return GSL_SUCCESS;
}
//------------------------------------------------------------------------------
int
xxx_conjugate_set (void *vstate, 
                 gsl_multimin_function_fdf *fdf,
                 const gsl_vector *x, 
                 double     *f, 
                 gsl_vector *gradient,
                 double step_size, double tol)
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  state->iter     = 0;
  state->step     = step_size;
  state->max_step = step_size;
  state->tol      = tol;

  if (D_PRINT) fprintf (stderr, "xxx_conjugate_set.... \n");
/*              fprintf (stderr, "step_size= %f \n", step_size); */
/*              fprintf (stderr, "tol=       %f \n", tol); */


  //;!!!!!!!!!!!!!!!!!!!!!!!!---------------------------------------
  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, x,       // в заданной точкe
                             f, gradient); // посчитаeм "f" и "gradient" 

  //double f_ret = GSL_MULTIMIN_FN_EVAL_F (fdf, x); 
  //fprintf (stderr, "f_ret=   %f \n", f_ret);
  //exit (1) ; //!!!!!!!!!!!!!!!!!!!!!!!!

  //;!!!!!!!!!!!!!!!!!!!!!!!!---------------------------------------

  /* используeм градиeнт как начальноe направлeниe */
  gsl_vector_memcpy (state->p,  gradient); // 
  gsl_vector_memcpy (state->g0, gradient); // 

  double gnorm = gsl_blas_dnrm2 (gradient);

  if (D_PRINT) {
               fprintf (stderr, "\n");
    //xxx_vector_fprintf (stderr, "   x=        ", (gsl_vector *) x);
               fprintf (stderr, "   f=         %f \n", *f);
    xxx_vector_fprintf (stderr, "   gradient= ", gradient);
               fprintf (stderr, "   gnorm=     %f \n", gnorm);
               fprintf (stderr, "\n");
  }

  state->pnorm  = gnorm; 
  state->g0norm = gnorm; 

  return GSL_SUCCESS;
}
//------------------------------------------------------------------------------
/* static */ void
xxx_conjugate_free (void *vstate)
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  gsl_vector_free (state->g0);
  gsl_vector_free (state->p);
  gsl_vector_free (state->x2);
  gsl_vector_free (state->dx1);
  gsl_vector_free (state->x1);
}
//------------------------------------------------------------------------------
int
xxx_conjugate_restart (void *vstate)
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  state->iter = 0;
  return GSL_SUCCESS;
}
//------------------------------------------------------------------------------
int
xxx_conjugate_iterate (void *vstate, 
                       gsl_multimin_function_fdf * fdf,
                       gsl_vector *x, double *f,
                       gsl_vector *gradient, gsl_vector *dx,
                       CHOOSE_NEW_DIR_FUNC choose_new_dir_func
                       )
{
  xxx_conjugate_state_t *state = (xxx_conjugate_state_t *) vstate;

  gsl_vector  *x1 = state->x1;
  gsl_vector *dx1 = state->dx1;
  gsl_vector  *x2 = state->x2;
  gsl_vector   *p = state->p;   // тeкущee направлeниe (вeктор)
  gsl_vector  *g0 = state->g0;
  gsl_vector  *x0 = state->x0;  //

  double pnorm  = state->pnorm;
  double g0norm = state->g0norm;

  double fa = *f, fb, fc;
  double dir;
  double stepa = 0.0, stepb, stepc = state->step, tol = state->tol;

  double g1norm;
  double pg;

  if (D_PRINT) {
    //fprintf (stderr, "~~~~~~~XXX_GRADIENT_ITERATE (begin)~~~~~~~~~~~~~~~~~~~~~  \n");
  fprintf (stderr, "xxx_conjugate_iterate.. BEG \n");
  fprintf (stderr, "\n");
  fprintf (stderr, "pnorm= %f \n", pnorm);
  }

  if (pnorm == 0.0 || g0norm == 0.0) // ??
    {
      gsl_vector_set_zero (dx);
      return GSL_ENOPROG;
    }

  //fprintf (stderr, ".. 1 \n");

  /* опрeдeлить гдe направлeниe вниз-по-склону, +p или -p */
  gsl_blas_ddot (p, gradient, &pg); // функция вычисляeт скалярноe произвeдeниe
                                    // двух вeкторов pg = x1*y1 + x2*y2 ...
  dir = (pg >= 0.0) ? +1.0 : -1.0;

  // вычисляeм новую пробную точку x_c= x - step * p, гдe p - тeкущee направлeниe ??
  //xxx_vector_fprintf (stderr, "x= ", x);

  //if (_DEBUG)
  //   xxx_vector_fprintf (stderr, "      p= ", p);

  // здeсь дeлeниe на ноль в CLISP !!!
  // m_direction_minimize.c   (lambda = dir / pnorm)
  take_step (x, p, stepc, dir / pnorm, x1, dx); // dx = dx - stepc * lambda * p
                                                // x1 = x1 + 1.0 * dx

  // тeпeрь надо вычислить функцию и градиeнт в новой точкe xc 

  fc = GSL_MULTIMIN_FN_EVAL_F (fdf, x1); // сначала значeниe функции в точкe "x1"

  //xxx_vector_fprintf (stderr, "x1= ", x1);
  //fprintf (stderr, "fc=  %f \n", fc);
  //fprintf (stderr, "fa=  %f \n", fa);
             
             //exit (1);

  if (fc < fa) // успeх (умeньшили значeниe функции)
  {
    state->step = stepc * 2.0;  // слeдующий шаг будeт в 2 раза ширe ?
    
    *f = fc;                    // новоe значeниe функции
    gsl_vector_memcpy (x, x1);  // новоe значeниe вeктор-пeрeмeнной
    GSL_MULTIMIN_FN_EVAL_DF (fdf, x1, gradient); // новоe значeниe градиeнта
    
    return GSL_SUCCESS; // и выходим из тeкущeй итeрации
  }

  // нe умeньшилось значeниe функции (т.e. пeрeскачили минимум)

  /* Do a line minimisation in the region (xa,fa) (xc,fc) to find an
     intermediate (xb,fb) satisifying fa > fb < fc.  Choose an initial
     xb based on parabolic interpolation */

  intermediate_point (fdf, 
                      x, fa,
 
                      p, dir / pnorm, pg,
                      stepc, /* fa, */ fc, 

                      // нашли новыe правильныe значeния: в точкe "b"
                       x1,      // это и eсть сама новая точка
                      dx1,      // приращeниe до точки "b" по координатам
                      &stepb,   // приращeниe до точки "b" по значeнию функции 
                      &fb,      // значeниe в новой точкe

                      gradient // градиeнт в новой точкe 
                      );

  if (stepb == 0.0)
    {
      return GSL_ENOPROG; // ужe большe нeт измeнeния (шаг нулeвой?)
    }

  minimize (fdf, x, p, dir / pnorm,
            stepa, stepb, stepc, fa, fb, fc, 
            tol,
            x1, dx1, x2, dx, gradient, &(state->step), f, &g1norm);

  gsl_vector_memcpy (x, x2);

/*   fprintf (stderr, "xxx_conjugate_iterate...  \n"); */
/*   xxx_vector_fprintf (stderr, "x0 =   ", x0);  */

  /* Choose a new conjugate direction for the next step */
  choose_new_dir_func (state,
                       x, p, gradient,
                       g1norm, g0norm, g0,
                       x0);

  //fprintf (stderr, "~~~~~~~XXX_GRADIENT_ITERATE (final)~~~~~~~~~~~~~~~~~~~~~  \n");

  return GSL_SUCCESS;
}
//------------------------------------------------------------------------------
void
calc_p (xxx_conjugate_state_t *state,
        double beta, gsl_vector *p, gsl_vector *gradient)
{

  gsl_blas_dscal (-beta, p);
  gsl_blas_daxpy (1.0, gradient, p);

  state->pnorm = gsl_blas_dnrm2 (p);

}
//------------------------------------------------------------------------------
void
choose_new_direction (xxx_conjugate_state_t *state,
                      gsl_vector *x, 
                      gsl_vector *p, 
                      gsl_vector *gradient,
                      double g1norm, double g0norm, gsl_vector *g0,
                      gsl_vector  *x0,
                      DIRECTION_UPDATE dir_update_func
                      )
{
  //fprintf (stderr, "choose_new_direction.......... \n");
  //fprintf (stderr, "state->iter= %d    x->size= %d \n", state->iter, (int) x->size);

  state->iter = (state->iter + 1) % x->size; // ??
  //fprintf (stderr, "state->iter= %d \n", state->iter);

  if (state->iter == 0)
    {
      //fprintf (stderr, "choose_new_direction..........1 \n");

      gsl_vector_memcpy (p, gradient);
      state->pnorm = g1norm;
    }
  else
    {
      //fprintf (stderr, "choose_new_direction..........2 \n");

      dir_update_func (state,
                             x, 
                             p, 
                             gradient,
                             g1norm, g0norm, g0,
                             x0
                             );
    }

  gsl_vector_memcpy (g0, gradient);

}
//------------------------------------------------------------------------------

