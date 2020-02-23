//------------------------------------------------------------------------------

double
calc_beta_pr (double g0norm, double g1norm,
              gsl_vector *gradient, gsl_vector *g0);

void
direction_update_pr (xxx_conjugate_state_t *state,
                       gsl_vector *x, 
                       gsl_vector *p, 
                       gsl_vector *gradient,
                       double g1norm, double g0norm, gsl_vector *g0,
                       gsl_vector  *x0
                     );

void
choose_new_direction_pr (xxx_conjugate_state_t *state,
                         gsl_vector *x, 
                         gsl_vector *p, 
                         gsl_vector *gradient,
                         double g1norm, double g0norm, gsl_vector *g0,
                         gsl_vector  *x0
                         );

//------------------------------------------------------------------------------
