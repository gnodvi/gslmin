
//------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------

extern int _DEBUG;
extern int D_PRINT; // для отладочной пeчати по BERGER


void
xxx_vector_fprintf_line (FILE *stream, char *name, gsl_vector *v);
void
xxx_vector_printf_line (char *name, gsl_vector *v);

void
xxx_vector_fprintf (FILE *stream, char *name, gsl_vector *v);
void
xxx_vector_printf (char *name, gsl_vector *v);

/* static */ void
take_step (const gsl_vector *x, const gsl_vector *p,
           double step, double lambda, 
           gsl_vector *x1, gsl_vector *dx);

//------------------------------------------------------------------------------

typedef struct {
  double step;
  double max_step;
  double tol;

  gsl_vector  *x1;
  gsl_vector  *g1;

  int    iter;

  gsl_vector *dx1;
  gsl_vector  *x2;

  double g0norm;
  double pnorm;

  gsl_vector   *p;
  gsl_vector  *g0;

  gsl_vector  *x0;      // vector_bfgs_state_t
  gsl_vector *dx0;
  gsl_vector *dg0;
}
xxx_conjugate_state_t;  


/* static */ int
xxx_conjugate_alloc (void *vstate, size_t n);

/* static */ int
xxx_conjugate_set (void *vstate, 
                 gsl_multimin_function_fdf *fdf,
                 const gsl_vector *x, 
                 double     *f, 
                 gsl_vector *gradient,
                 double step_size, double tol);

/* static */ void
xxx_conjugate_free (void *vstate);

int
xxx_conjugate_restart (void *vstate);

typedef void (*CHOOSE_NEW_DIR_FUNC) (xxx_conjugate_state_t *state,
                           gsl_vector *x, gsl_vector *p, 
                           gsl_vector *gradient,
                           double g1norm, double g0norm, gsl_vector *g0,
                           gsl_vector  *x0
                                  );

/* static */ int
xxx_conjugate_iterate (void *vstate, 
                       gsl_multimin_function_fdf * fdf,
                       gsl_vector *x, double *f,
                       gsl_vector *gradient, gsl_vector *dx,
                       CHOOSE_NEW_DIR_FUNC choose_new_dir_func
                       );

void
calc_p (xxx_conjugate_state_t *state,
        double beta, gsl_vector *p, gsl_vector *gradient);

typedef void
(*DIRECTION_UPDATE) (xxx_conjugate_state_t *state,
                       gsl_vector *x, 
                       gsl_vector *p, 
                       gsl_vector *gradient,
                       double g1norm, double g0norm, gsl_vector *g0,
                       gsl_vector  *x0
                     );

void
choose_new_direction (xxx_conjugate_state_t *state,
                      gsl_vector *x, 
                      gsl_vector *p, 
                      gsl_vector *gradient,
                      double g1norm, double g0norm, gsl_vector *g0,
                      gsl_vector  *x0,
                      DIRECTION_UPDATE dir_update_func
                      );

//------------------------------------------------------------------------------
#ifdef __cplusplus
}        
#endif
//------------------------------------------------------------------------------
