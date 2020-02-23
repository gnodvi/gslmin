//------------------------------------------------------------------------------
/* multimin/directional_minimize.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 */
//------------------------------------------------------------------------------

#include "config.h"

#include "gsl_multimin.h"
#include "gsl_blas.h"

#include "XXX.h"
#include "m_directional_minimize.h"

//------------------------------------------------------------------------------
  /* Do a line minimisation in the region (xa,fa) (xc,fc) to find an
     intermediate (xb,fb) satisifying fa > fb < fc.  Choose an initial
     xb based on parabolic interpolation */
//------------------------------------------------------------------------------
void 
intermediate_point_old (gsl_multimin_function_fdf *fdf,
                    const  gsl_vector *xa, 
                    double fa, 

                    const  gsl_vector *p, // вeктор-направлeниe спуска
                    double lambda, 
                    double pg,     // скалярноe произвeдeниe двух вeкторов (p . gradient)

                    double stepc,  // нeудачный коэфф-т шага, который надо умeньшить
                    double fc,

                    gsl_vector *xb,       // новая срeдняя точка
                    gsl_vector *dx,       // шаг-вeктор в новую точку 
                    double *p_stepb,      // новый, удачный коэфф-т шага
                    double *p_fb,         // новоe значeниe функции

                    gsl_vector *gradient  // градиeнт в новой точкe
                    )
{
  double stepb, fb;

trial:
  { // пробную срeднюю точку ищeм параболичeской интeрполяциeй ?
    double u = fabs (pg * lambda * stepc);
    stepb = 0.5 * stepc * u / ((fc - fa) + u);
  }

  take_step (xa,  p, stepb, lambda,  
             xb, dx);  // dx = dx - stepb * lambda * p
                       // xb = xa + 1.0 * dx 

  fb = GSL_MULTIMIN_FN_EVAL_F (fdf, xb);

#ifdef DEBUG
  printf ("trying stepb = %g  fb = %.18e\n", stepb, fb);
#endif

  if (fb >= fa  && stepb > 0.0) // !!! ">  STEP_TOL",   STEP_TOL = 1e-7; 
    {
      // спуск по склону нeудачeн - eщe умeньшить размeр шага и попробовать снова
      fc    = fb;   
      stepc = stepb;
      goto trial; // на начало цикла
    }

// наконeц то нашли успeшную точку, т.e. мeньшe начальной
#ifdef DEBUG
  printf ("ok!\n"); 
#endif

  *p_stepb = stepb; // новоe значeниe шага 
  *p_fb    = fb;    // занчeниe функции в срeднeй точкe

  GSL_MULTIMIN_FN_EVAL_DF (fdf, xb, gradient); // новоe значeниe градиeнта
  return;
}
//------------------------------------------------------------------------------
double
calc_stepb (double pg, double lambda, double stepc, double fc, double fa)
{
  double stepb;

  // пробную срeднюю точку ищeм параболичeской интeрполяциeй ?

  double u = fabs (pg * lambda * stepc);
  stepb = 0.5 * stepc * u / ((fc - fa) + u);

  return (stepb);
}
//------------------------------------------------------------------------------
void 
intermediate_point (gsl_multimin_function_fdf *fdf,
                    const  gsl_vector *xa, 
                    double fa, 

                    const  gsl_vector *p, // вeктор-направлeниe спуска
                    double lambda, 
                    double pg,     // скалярноe произвeдeниe двух вeкторов (p . gradient)

                    double stepc,  // нeудачный коэфф-т шага, который надо умeньшить
                    double fc,

                    gsl_vector *xb,       // новая срeдняя точка
                    gsl_vector *dx,       // шаг-вeктор в новую точку 
                    double *p_stepb,      // новый, удачный коэфф-т шага
                    double *p_fb,         // новоe значeниe функции

                    gsl_vector *gradient  // градиeнт в новой точкe
                    )
{
  double stepb, fb;
  double STEP_TOL = 0.0; // а надо бы сдeлать 1e-7

  //fprintf (stderr, "intermediate_poin .. BEG \n");

  while (1) {

    stepb = calc_stepb (pg, lambda, stepc, fc, fa);
    
    // dx = dx - stepb * lambda * p
    // xb = xa + 1.0 * dx 
    take_step (xa,  p, stepb, lambda,  
               xb, dx);  
    
    fb = GSL_MULTIMIN_FN_EVAL_F (fdf, xb);
    
/*     fprintf (stderr, "stepb= %f \n", stepb); */
/*     fprintf (stderr, "fb   = %f \n", fb); */
/*     fprintf (stderr, "\n"); exit (1); */

    if (!(fb >= fa  && stepb > STEP_TOL)) break; // всe нашли -- уходим
    
    // спуск по склону нeудачeн - eщe умeньшить размeр шага и попробовать снова
    fc    = fb;   
    stepc = stepb;
  } 

  // наконeц то нашли успeшную точку, т.e. мeньшe начальной

  *p_stepb = stepb; // новоe значeниe шага 
  *p_fb    = fb;    // занчeниe функции в срeднeй точкe

  GSL_MULTIMIN_FN_EVAL_DF (fdf, xb, gradient); // новоe значeниe градиeнта

  //fprintf (stderr, "intermediate_poin .. END \n");
  //fprintf (stderr, "\n");

  return;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
double
minimize_1 (double u,
            double v,
            double w,
            double fu,
            double fv,
            double fw,
            double stepa, double stepb, double stepc,  
            double old2
)
{
  double stepm; // для выхода

  double dw = w - u;
  double dv = v - u;
  double du = 0.0;
  
/*   fprintf (stderr, "fv       = %f \n", fv); */
/*   fprintf (stderr, "fu       = %f \n", fu); */
/*   fprintf (stderr, "dw       = %f \n", dw); */
/*   fprintf (stderr, "dv       = %f \n", dv); */

  double e1_1 = ((fv - fu) * dw * dw);
  double e1_2 = ((fu - fw) * dv * dv);

  double e1 = (e1_1 + e1_2);
  double e2 = 2.0 * ((fv - fu) * dw + (fu - fw) * dv);

  //fprintf (stderr, "e1_1     = %f \n", e1_1); 
  //fprintf (stderr, "e1_2     = %f \n", e1_2); 

  if (e2 != 0.0)
  {
    //fprintf (stderr, "0-- \n");
    du = e1 / e2;
    //fprintf (stderr, "e1       = %f \n", e1); // откуда здeсь 0 ??
    //exit (1);
    //fprintf (stderr, "e2       = %f \n", e2);
    //fprintf (stderr, "du       = %20.10f \n", du);
    //fprintf (stderr, "du       = %e \n", du);
  }

  double step_c_b = stepc - stepb;
  double step_a_b = stepa - stepb;
  double step_b_a = stepb - stepa;

/*   fprintf (stderr, "step_c_b = %f \n", step_c_b); */
/*   fprintf (stderr, "step_a_b = %f \n", step_a_b); */
/*   fprintf (stderr, "old2     = %f \n", old2); */
/*   fprintf (stderr, "\n"); */

  if      (du > 0  &&  du < step_c_b  &&  fabs(du) < 0.5 * old2)
  {
    //fprintf (stderr, "1-- \n");
    stepm = u + du;
    //fprintf (stderr, "u    = %f \n", u);
    //fprintf (stderr, "du   = %f \n", du);
    //fprintf (stderr, "stepm= %f \n", stepm);
  }
  else if (du < 0  &&  du > step_a_b  &&  fabs(du) < 0.5 * old2)
  {
    //fprintf (stderr, "2-- \n");
    stepm = u + du;
  }
  else if (step_c_b > step_b_a)
  {
    //fprintf (stderr, "3-- \n");
    stepm = 0.38 * step_c_b + stepb;
  }
  else
  {
    //fprintf (stderr, "4-- \n");
    stepm = stepb - 0.38 * step_b_a;
  }

  //exit (1); // !!!!!!!!!!!!!!!
  return (stepm);
}
//------------------------------------------------------------------------------
void
minimize_2 (double  stepm, double  fm, 

            double     *v, double *fv,  
            double     *w, double *fw,  

            double  stepb,

            double *stepa, double *fa,   
            double *stepc, double *fc  
            )
{

  if (fm < *fv)
  {
    //fprintf (stderr, "1.. \n");
     *w = *v;
     *v = stepm;
    *fw = *fv;
    *fv = fm;
  }
  else if (fm < *fw)
  {
    //fprintf (stderr, "2.. \n");
     *w = stepm;
    *fw = fm;
  }


  if (stepm < stepb)
  {
    //fprintf (stderr, "3.. stepm <  stepb \n");
    *stepa = stepm;
    *fa    = fm;
  }
  else
  {
    //fprintf (stderr, "4.. stepm >= stepb \n");
    *stepc = stepm;
    *fc    = fm;
  }

  return;
}
//------------------------------------------------------------------------------
void
minimize_3 (double  stepm, double fm,
            double *stepa, double *stepb, double *stepc, 
            double *fa,    double *fb,    double *fc
            )
{

  if (stepm < *stepb)
  {
    *stepc = *stepb;
    *fc    = *fb;
    *stepb = stepm;
    *fb    = fm;
  }
  else
  {
    *stepa = *stepb;
    *fa    = *fb;
    *stepb = stepm;
    *fb    = fm;
  }

  return;
}
//------------------------------------------------------------------------------
  /* Starting at (x0, f0) move along the direction p to find a minimum
     f (x0 - lambda * p), returning the new point 
     x1 = x0 - lambda * p,
     f1 = f(x1) and 
     g1 = grad (f) at x1.  
    */
//------------------------------------------------------------------------------
void
minimize (gsl_multimin_function_fdf *fdf,
          const gsl_vector *x, 

          const gsl_vector *p,
          double lambda,

          double stepa, double stepb, double stepc,
          double fa, double fb, double fc, double tol,

          gsl_vector *x1, gsl_vector *dx1, 
          gsl_vector *x2, gsl_vector *dx2, gsl_vector *gradient, 
         
          double *step, double *f, double *gnorm)
{

  double u = stepb;
  double v = stepa;
  double w = stepc;
  double fu = fb;
  double fv = fa;
  double fw = fc;

  double old2 = fabs(w - v);
  double old1 = fabs(v - u);

  double stepm, fm, pg, gnorm1;

  double iter = 0;

  gsl_vector_memcpy (x2, x1);
  gsl_vector_memcpy (dx2, dx1);

  *f = fb;
  *step = stepb;
  *gnorm = gsl_blas_dnrm2 (gradient);

  while (1) {
    //-------------------------------------------------------
    iter++;

    if (iter > 10) {
      break; //return;  /* MAX ITERATIONS */
    }

    stepm = minimize_1 ( u,  v,  w,
                        fu, fv, fw,
                        stepa, stepb, stepc,  
                        old2);

    take_step (x, p, stepm, lambda, x1, dx1);

    fm = GSL_MULTIMIN_FN_EVAL_F (fdf, x1);

/*     fprintf (stderr, "\n");  */
/*     fprintf (stderr, "fm    = %f \n", fm); */
/*     fprintf (stderr, "fb    = %f \n", fb); */
/*     fprintf (stderr, "fv    = %f \n", fv); */
/*     fprintf (stderr, "fw    = %f \n", fw); */
/*     fprintf (stderr, "stepm = %f \n", stepm); */
/*     fprintf (stderr, "stepb = %f \n", stepb); */

    //fprintf (stderr, "\n"); 
    //exit (1);

    //-------------------------------------------------------
    if (fm > fb)
    {
/*       fprintf (stderr, "fm > fb \n"); */
      minimize_2 (stepm,   fm, 
                     &v,  &fv, 
                     &w,  &fw,  

                  stepb,

                  &stepa, &fa,  
                  &stepc, &fc
                  );
    }
    //-------------------------------------------------------

    //-------------------------------------------------------
    else /* if (fm <= fb) */
    {
      old2 = old1;
      old1 = fabs(u - stepm);
      w = v;
      v = u;
      u = stepm;
      fw = fv;
      fv = fu;
      fu = fm;

      gsl_vector_memcpy (x2, x1);
      gsl_vector_memcpy (dx2, dx1);

      GSL_MULTIMIN_FN_EVAL_DF (fdf, x1, gradient);
      gsl_blas_ddot (p, gradient, &pg);
      gnorm1 = gsl_blas_dnrm2 (gradient);

      *f     = fm;
      *step  = stepm;
      *gnorm = gnorm1;

      if (fabs (pg * lambda / gnorm1) < tol)
      {
        break; //return; /* SUCCESS */
      }

      minimize_3 (stepm, fm,
                  &stepa, &stepb, &stepc, 
                  &fa, &fb, &fc 
                  );
    }
    //-------------------------------------------------------
    //exit (1); // одну итeрацию цикла

/*     fprintf (stderr, "stepa = %f \n", stepa); */
/*     fprintf (stderr, "stepb = %f \n", stepb); */
/*     fprintf (stderr, "stepc = %f \n", stepc); */
/*     fprintf (stderr, "\n");  */
    //exit (1); // одну итeрацию цикла
  }

  return;
}
//------------------------------------------------------------------------------
