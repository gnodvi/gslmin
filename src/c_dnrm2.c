// -*-  mode: c    ; coding: koi8   -*- ----------------------------------------

//------------------------------------------------------------------------------

#include "gsl_math.h"
#include "gsl_cblas.h"

#include "cblas.h"

double
cblas_dnrm2 (const int N, const double *X, const int incX)
{
#define BASE double
#include "c_source_nrm2_r.h"
#undef BASE
}
