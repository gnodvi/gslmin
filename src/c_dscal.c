// -*-  mode: c    ; coding: koi8   -*- ----------------------------------------

//------------------------------------------------------------------------------

#include "gsl_math.h"
#include "gsl_cblas.h"

#include "cblas.h"

void
cblas_dscal (const int N, const double alpha, double *X, const int incX)
{
#define BASE double
#include "c_source_scal_r.h"
#undef BASE
}
