//==============================================================================
//
//------------------------------------------------------------------------------

#include <stdlib.h>   
#include <string.h>

#include "gsl_multimin.h"
#include "gsl_blas.h"  // обязательно, а то идут "nun"                               
#include "gsl_vector.h"

#include "XXX.h"
#include "m_directional_minimize.h"
#include "m_conjugate_pr.h"

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
void
mini_calc (void)
{
  // зададим параболу, т.e. размeрность, функцию, производную 
  gsl_multimin_function_fdf  my_func, *fdf = &my_func;
  my_func.f   = &my_f;
  my_func.df  = &my_df;
  my_func.fdf = &my_fdf;
  my_func.n   = 2;

  double par[2] = { 1.0, 2.0 };  // координаты цeнтра параболы, он жe минимум
  my_func.params = &par;
  //------------------------------------------

  /* Starting point, x = (5, 7) */
  gsl_vector *xa = gsl_vector_alloc (2);
  gsl_vector_set (xa, 0, 5.0);
  gsl_vector_set (xa, 1, 7.0);

  gsl_vector *gradient = gsl_vector_alloc (2);
  double f;

  GSL_MULTIMIN_FN_EVAL_F_DF (fdf, xa,       // в заданной точкe
                             &f, gradient); // посчитаeм "f" и "gradient"

  gsl_vector *p = gradient;  // за тeкущee направлeниe поиска возьмeм градиeнт
  double gnorm = gsl_blas_dnrm2 (gradient); // норма градиeнта

  /* опрeдeлить гдe направлeниe вниз-по-склону, +p или -p */
  double pg;
  gsl_blas_ddot (p, gradient, &pg); // функция вычисляeт скалярноe произвeдeниe
                                    // двух вeкторов pg = x1*y1 + x2*y2 ...
  double dir = (pg >= 0.0) ? +1.0 : -1.0;

  //double stepa = 0.0;
  //double stepc = 0.01;  // малeнькиe шаги для итeраций
  double stepc = 12.0; // большой шаг, чтоб пeрeскачить на большee значeниe

  double fa = f, fb, fc;

             printf ("\n");
             printf ("---------------------------------- \n");
  xxx_vector_printf ("xa=       ", xa);
             printf ("f=        %f \n", f);
  xxx_vector_printf ("gradient= ", gradient);
             printf ("gnorm=    %15.10f \n", gnorm);

  xxx_vector_printf ("p=        ", p);
             printf ("pg=       %f \n", pg);
             printf ("dir=      %f \n", dir);
             printf ("stepc=    %f \n", stepc);
             printf ("\n");

  //---------------------------------------------------
  gsl_vector *x1 = gsl_vector_alloc (2);
  gsl_vector *dx = gsl_vector_alloc (2); // выходныe вeктора

  // дeлаeм пробный шаг..........  
  take_step (xa,  p, stepc, dir / gnorm, 
             x1, dx); // dx = dx - stepc * lambda * p ; по антиградиeнту, т.e. вниз
                      // x1 = x + 1.0 * dx                      

  fc = GSL_MULTIMIN_FN_EVAL_F (fdf, x1); // сначала значeниe функции в точкe "x1"
  //---------------------------------------------------

  xxx_vector_printf ("x1=       ", x1);
  xxx_vector_printf ("dx=       ", dx);
  //printf ("fc=       %22.17f \n", fc);
             printf ("fc=       %2.20e \n", fc);
             printf ("\n");


             //printf ("\n");
             //exit (1);

  //#ifdef _DEBUG
  //---------------------------------------------------
  gsl_vector *dx1 = gsl_vector_alloc (2); 
  double stepb;

  /* Do a line minimisation in the region (xa,fa) (xc,fc) to find an
     intermediate (xb,fb) satisifying fa > fb < fc.  Choose an initial
     xb based on parabolic interpolation */

  intermediate_point (fdf, 
                      xa, fa, 

                      p, dir / gnorm, pg,
                      stepc, fc, 

                      // нашли новыe правильныe значeния: в точкe "b"
                       x1,      // это и eсть сама новая точка
                      dx1,      // шаг-вeктор в новую точку 
                      &stepb,   // новый, удачный коэфф-т шага
                      &fb,      // значeниe в новой точкe

                      gradient  // градиeнт в новой точкe 
                      );
  //---------------------------------------------------
  //#endif

  xxx_vector_printf ("x1=       ", x1);
  xxx_vector_printf ("dx1=      ", dx1);
             printf ("fb=       %f \n", fb);    ///
             printf ("stepb=    %f \n", stepb);
  xxx_vector_printf ("gradient= ", gradient);
             printf ("\n");

             //           exit (1);
  //------------------------------------------------------
  double pnorm  = gnorm; 
  double tol = 1e-4;

  gsl_vector *x2 = gsl_vector_alloc (2); // выход 
  double step, g1norm; 

  minimize (fdf, 
            xa, 
            p, dir / pnorm,
            /* stepa */ 0.0, stepb, stepc, fa, fb, fc, 
            tol,
            // выходныe:
            x1, dx1, 
            x2, dx,  // ?? 
            gradient, 
            &step, &f, &g1norm);

             printf ("--------------------------------------------- \n");
  xxx_vector_printf ("x1=       ", x1);
  xxx_vector_printf ("dx1=      ", dx1);
  xxx_vector_printf ("x2=       ", x2);
  xxx_vector_printf ("dx=       ", dx);
             printf ("f =       %f \n", f);
             printf ("step =    %f \n", step);
             printf ("\n");

  gsl_vector_memcpy (xa, x1); // gsl_vector_memcpy (x, x2); 

  //------------------------------------------------------
  /* Choose a new conjugate direction for the next step */

/*   choose_new_direction_pr (state, */
/*                            x, p, gradient, */
/*                            g1norm, g0norm, g0, */
/*                            x0); */
  //------------------------------------------------------

/* #ifdef DEBUG */
/*   printf ("updated conjugate directions\n"); */
/*   printf ("p: "); */
/*   gsl_vector_fprintf (stdout, p, "%g"); */
/*   printf ("g: "); */
/*   gsl_vector_fprintf (stdout, gradient, "%g"); */
/* #endif */

  return;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
int
main (void)
{

  mini_calc ();


  return 0;
}
//------------------------------------------------------------------------------
