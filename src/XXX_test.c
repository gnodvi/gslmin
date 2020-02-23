//------------------------------------------------------------------------------

//#include <config.h>

#include <stdlib.h>

#include "gsl_test.h"
#include "gsl_blas.h"
#include "gsl_multimin.h"

//#include <gsl/gsl_ieee_utils.h>
//#include "Test_funcs.h"

#include "XXX.h"

//------------------------------------------------------------------------------
int
main (void)
{

  gsl_vector *x = gsl_vector_alloc (2);

  gsl_vector_set (x, 0, 5.0);
  gsl_vector_set (x, 1, 7.0);

  xxx_vector_fprintf (stderr, "x= ", x);

  exit (0);
}
//------------------------------------------------------------------------------
