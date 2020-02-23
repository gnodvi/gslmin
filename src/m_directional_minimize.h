//------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------

/* static */ void 
intermediate_point (gsl_multimin_function_fdf * fdf,
                    const  gsl_vector *x,
                    double fa,
                    
                    const gsl_vector *p,
                    double lambda, 
                    double pg,     // ��������e ������e�e��e ���� �e������ (p . gradient)

                    double stepc,  // �e������� ���, ������� ���� ��e������
                    /* double fa, */ double fc,

                    gsl_vector *x1,       // ����� �����
                    gsl_vector *dx,       // ���-�e���� ������e���
                    double *step,         // �������������� ��� � ����� �����
                    double *f,            // ����e ����e��e �������

                    gsl_vector *gradient // �����e�� � ����� ����e
                    );

/* static */ void
minimize (gsl_multimin_function_fdf * fdf,
          const gsl_vector * x, const gsl_vector * p,
          double lambda,
          double stepa, double stepb, double stepc,
          double fa, double fb, double fc, double tol,
          gsl_vector * x1, gsl_vector * dx1, 
          gsl_vector * x2, gsl_vector * dx2, gsl_vector * gradient,          
          double * step, double * f, double * gnorm);

//------------------------------------------------------------------------------
#ifdef __cplusplus
}        
#endif
//------------------------------------------------------------------------------
