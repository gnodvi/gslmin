// -*-  mode: c    ; coding: koi8   -*- ----------------------------------------

//------------------------------------------------------------------------------

#include "gsl_math.h"
#include "gsl_cblas.h"

#include "cblas.h"

void
cblas_daxpy (const int N, const double alpha, const double *X, const int incX,
             double *Y, const int incY)
{
#define BASE double
#include "c_source_axpy_r.h"
#undef BASE
}
