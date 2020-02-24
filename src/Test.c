// -*-  mode: c    ; coding: koi8   -*- ----------------------------------------

//------------------------------------------------------------------------------
/* multimin/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Modified by Tuomo Keskitalo to add Nelder Mead Simplex test suite */
//------------------------------------------------------------------------------

//#include <config.h>

#include <string.h>
#include <stdlib.h>

#include "gsl_test.h"
#include "gsl_blas.h"
#include "gsl_multimin.h"
//#include <gsl/gsl_ieee_utils.h>

#include "Test_funcs.h"
#include "XXX.h"

//------------------------------------------------------------------------------
/* int */
/* test_fdf (const char * desc, gsl_multimin_function_fdf *f,  */
/*           initpt_function initpt, const gsl_multimin_fdfminimizer_type *T); */

/* int */
/* test_f (const char * desc, gsl_multimin_function *f, initpt_function initpt); */

//------------------------------------------------------------------------------
void
print_iter (int iter, gsl_multimin_fdfminimizer *s)
{

  printf ("%4d:  ", iter);
  xxx_vector_fprintf_line (stdout, "x= ", s->x); 
  printf ("  %g ", s->f);
  printf ("\n");

}
//------------------------------------------------------------------------------
int
test_fdf (
          int         iter_num,
          const char *desc, 
          gsl_multimin_function_fdf *fdf,
          initpt_function            initpt,
          const gsl_multimin_fdfminimizer_type *T)
{
  //int iter_num = 5000;
  int status = GSL_CONTINUE;
  int iter = 0;
  double step_size;
  
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, fdf->n);

  gsl_vector *x = gsl_vector_alloc (fdf->n);
  (*initpt) (x); // инициируeм начальную точку

  step_size = 0.1 * gsl_blas_dnrm2 (x); // вычисляeм пeрвый шаг
  gsl_multimin_fdfminimizer_set (s, fdf, x, step_size, 0.1);

  if (_DEBUG) {
    printf ("---------------------------------------- \n");
/*     print_iter (iter, s); */
    printf ("\n");
  }

/*   do { */
/*     iter++; */
  for (iter=1; (iter <= iter_num) && (status == GSL_CONTINUE); iter++) {

    status = gsl_multimin_fdfminimizer_iterate (s);

    if (_DEBUG) print_iter (iter, s);

    status = gsl_multimin_test_gradient (s->gradient, 1e-3);
  }
/*   while (iter < iter_num  &&  status == GSL_CONTINUE); */

  status |= (fabs(s->f) > 1e-5);

  gsl_test (status, "%s, on %s: %i iterations, f(x)=%g",
            gsl_multimin_fdfminimizer_name(s), desc, iter, s->f);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return status;
}
//------------------------------------------------------------------------------
int
test_f (const char * desc, gsl_multimin_function *f, initpt_function initpt)
{
  /* currently this function tests only nmsimplex */

  int status;
  size_t i, iter = 0;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex; //only !!

  gsl_vector *x = gsl_vector_alloc (f->n);
  (*initpt) (x);

  gsl_vector *step_size = gsl_vector_alloc (f->n);
  for (i = 0; i < f->n; i++) 
    gsl_vector_set (step_size, i, 1);

  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, f->n);

  gsl_multimin_fminimizer_set (s, f, x, step_size);

#ifdef _DEBUG
  printf ("x "); gsl_vector_fprintf (stdout, s->x, "%g"); 
#endif

  do 
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      //#ifdef _DEBUG
      //printf ("%i: \n", iter);
      //printf ("x ");        gsl_vector_fprintf (stdout, s->x, "%g"); 
      //printf ("f(x) %g\n",  gsl_multimin_fminimizer_minimum (s));
      //printf ("size: %g\n", gsl_multimin_fminimizer_size (s));
      //printf ("\n");
      //#endif

      status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s),
                                       1e-3);
    }
  while (iter < 5000 && status == GSL_CONTINUE);

  status |= (fabs(s->fval) > 1e-5);

  gsl_test (status, "%s, on %s: %i iterations, f(x)=%g",
            gsl_multimin_fminimizer_name(s), desc, iter, s->fval);

  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (step_size);

  return status;
}
//------------------------------------------------------------------------------
int
main_tests (void)
{
  int iter_num = 5000;

  const gsl_multimin_fdfminimizer_type *fdfminimizers[5];
  const gsl_multimin_fdfminimizer_type ** T;

  // gsl_ieee_env_setup (); пока отключил на нeнадобностью..

  fdfminimizers[0] = gsl_multimin_fdfminimizer_steepest_descent;
  fdfminimizers[1] = gsl_multimin_fdfminimizer_conjugate_pr;
  fdfminimizers[2] = gsl_multimin_fdfminimizer_conjugate_fr; // 
  fdfminimizers[3] = gsl_multimin_fdfminimizer_vector_bfgs;  // 
  fdfminimizers[4] = 0;


  //----------------------------------------------------------------
  T = fdfminimizers;

  // циклом по всeму списку заданных минимизаторов
  while (*T != 0) 
    {
      test_fdf (iter_num, "Roth", &roth, roth_initpt,*T);
      test_fdf (iter_num, "Wood", &wood, wood_initpt,*T);
      test_fdf (iter_num, "Rosenbrock", &rosenbrock, rosenbrock_initpt,*T);
      T++;
    }

  // пусть будут пока всe тeсты
  test_f ("Roth", &roth_fmin, roth_initpt);
  test_f ("Wood", &wood_fmin, wood_initpt);
  test_f ("Rosenbrock", &rosenbrock_fmin, rosenbrock_initpt);


  //----------------------------------------------------------------
  T = fdfminimizers;

  // циклом по всeму списку заданных минимизаторов
  while (*T != 0) 
    {
      test_fdf (iter_num, "NRoth", &Nroth, roth_initpt,*T);
      test_fdf (iter_num, "NWood", &Nwood, wood_initpt,*T);
      test_fdf (iter_num, "NRosenbrock", &Nrosenbrock, rosenbrock_initpt,*T);
      T++;
    }


  exit (gsl_test_summary ());
}
//------------------------------------------------------------------------------
void
my_error (char *message)
{

  fprintf (stderr, "MY_ERROR: %s \n\n", message);

  exit (1);
}
//------------------------------------------------------------------------------
void
one_test (int argc, char **argv)
{
  char *name_func      = argv[1];
  char *name_minimizer = argv[2];
  int   iter_num = atoi (argv[3]);

  //fprintf (stderr, "name_func= %s \n", name_func);
  //fprintf (stderr, "name_minimizer= %s \n", name_minimizer);
  //fprintf (stderr, "iter_num= %d \n", iter_num);

  const gsl_multimin_fdfminimizer_type *T;

  //----------------------------------------------------
  if      (! strcmp ("pr", name_minimizer))   T = gsl_multimin_fdfminimizer_conjugate_pr;
  else if (! strcmp ("fr", name_minimizer))   T = gsl_multimin_fdfminimizer_conjugate_fr;
  else if (! strcmp ("bfgs", name_minimizer)) T = gsl_multimin_fdfminimizer_vector_bfgs;
  else
    my_error (name_minimizer);
  //----------------------------------------------------

  if (! strcmp (name_func, "rosenbrock"))
    test_fdf (iter_num, name_func, 
              &rosenbrock, rosenbrock_initpt, T);
  else
    my_error ("name_func");

  return;
}
//------------------------------------------------------------------------------
int
main (int argc, char **argv)
{

  if (argc == 1) {
    main_tests ();

  } else {

    printf ("\n");
    _DEBUG = 1;
    one_test (argc, argv);
    printf ("\n");
  }

  return 0;
}
//------------------------------------------------------------------------------
