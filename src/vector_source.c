/* vector/vector_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 */

#ifndef HIDE_INLINE_STATIC

//------------------------------------------------------------------------------
BASE
FUNCTION (gsl_vector, get) (const TYPE (gsl_vector) * v, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          BASE zero = ZERO;
          GSL_ERROR_VAL ("1: index out of range", GSL_EINVAL, zero);
        }
    }

  /* The following line is a generalization of return v->data[i] */

  return *(BASE *) (v->data + MULTIPLICITY * i * v->stride);
}
//------------------------------------------------------------------------------
//void
//FUNCTION (gsl_vector, set) (TYPE (gsl_vector) * v, const size_t i, BASE x)
//------------------------------------------------------------------------------
void
gsl_vector_set (gsl_vector *v, const size_t i, double x)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          GSL_ERROR_VOID ("2: index out of range", GSL_EINVAL);
        }
    }

  //fprintf (stderr, "c_gsl_vector_set: i= %ld x= %f \n", i, x);

  /* The following line is a generalization of v->data[i] = x */

  *(double *) (v->data + MULTIPLICITY * i * v->stride) = x;
  //  *(BASE *) (v->data + MULTIPLICITY * i * v->stride) = x;
}
//------------------------------------------------------------------------------
BASE *
FUNCTION (gsl_vector, ptr) (TYPE (gsl_vector) * v, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          GSL_ERROR_NULL ("3: index out of range", GSL_EINVAL);
        }
    }

  return (BASE *) (v->data + MULTIPLICITY * i * v->stride);
}
//------------------------------------------------------------------------------
const BASE *
FUNCTION (gsl_vector, const_ptr) (const TYPE (gsl_vector) * v, const size_t i)
{
  if (gsl_check_range)
    {
      if (i >= v->size)         /* size_t is unsigned, can't be negative */
        {
          GSL_ERROR_NULL ("4: index out of range", GSL_EINVAL);
        }
    }

  /* The following line is a generalization of return v->data[i] */

  return (const BASE *) (v->data + MULTIPLICITY * i * v->stride);
}
//------------------------------------------------------------------------------

#endif
