# -*-   mode: makefile ; coding: koi8   -*- ------------------------------------

# =============================================================================# 

#GCC_o64   = gcc -D_OPT64
#GCC_p32   = gcc -malign-double -Wall
#CC    = gcc -Wall -g -static $(GLIB_INCS)

CC    = gcc -Wall -fpic # -g $(GLIB_INCS)

#CPP   = c++ -O2
#FLAGS = $(OPTIM) -I. -I./Others -D_UNIX 
#FLAGS = -D_UNIX 

GSL_SRC = /home/.../gsl_new/1.8/.../gsl-1.8
GSL_DIR = /home/.../gsl_new/1.8/.../R

MINGSL_SO = libgslmin.so  # 
MYGSL_SO  = libmygsl.so   # оболочка для вызова чeрeз нe GSL (libgslmin.so) из Лиспа 

#MINGSL_SO = libgslmin.dll  # 
#MYGSL_SO  = libmygsl.dll   # оболочка для вызова чeрeз нe GSL (libgslmin.so) из Лиспа 

MINGSL_OBJS = L/m_fdfminimizer.o L/m_conjugate_fr.o L/m_convergence.o L/m_vector_bfgs.o L/m_fminimizer.o L/m_simplex.o L/m_steepest_descent.o L/m_conjugate_pr.o L/m_diff.o L/m_directional_minimize.o \
	L/a_error.o L/a_stream.o L/a_results.o \
	L/vector.o L/v_init.o L/v_copy.o L/v_minmax.o \
	L/x_matrix.o L/x_init.o  L/x_getset.o \
	L/block.o  L/b_init.o    L/XXX.o \
	L/c_blas.o L/c_ddot.o L/c_dnrm2.o L/c_daxpy.o L/c_dscal.o 
#	\
#	L/MyGSL.o

# =============================================================================#

PROGS = Test Mini XXX_test Prim $(MINGSL_SO) #$(MYGSL_SO) 

all : $(PROGS) 

clean :
	rm -f  *.so *.dll L/*.o T/* $(PROGS) #T/*.o T/*.a 

# =============================================================================#

$(MINGSL_SO) : $(MINGSL_OBJS)
	gcc -shared -fpic -o $@ $^  

$(MYGSL_SO) : L/MyGSL.o
	gcc -shared -fpic -o $@ $^  

L/MyGSL.o : MyGSL.c 
	$(CC) -o $@ -c  $< #$(GSL_INC) 

# =============================================================================#

Test : Test.c gsl_test.h gsl_blas.h gsl_multimin.h Test_funcs.h T/test_funcs.o $(MINGSL_OBJS)
	$(CC) -o $@  $< T/test_funcs.o $(MINGSL_OBJS) -lm

T/test_funcs.o : Test_funcs.c Test_funcs.h gsl_multimin.h
	$(CC)  -o $@ -c  $< #$(GSL_INC) 

Prim : Prim.c  gsl_multimin.h gsl_blas.h gsl_vector.h $(MINGSL_SO)
	$(CC) -o $@  $< -L. -lgslmin -lm

Mini : Mini.c  gsl_multimin.h gsl_blas.h gsl_vector.h $(MINGSL_SO)
	$(CC) -o $@  $< -L. -lgslmin -lm

# =============================================================================#

GRADIENT_H = config.h gsl_multimin.h gsl_blas.h XXX.h m_directional_minimize.h

XXX_test : XXX_test.c  gsl_test.h gsl_blas.h gsl_multimin.h XXX.h $(MINGSL_OBJS)
	$(CC) -o $@  $< $(MINGSL_OBJS) -lm

L/XXX.o : XXX.c $(GRADIENT_H)
	$(CC) -o $@ -c  $< 

# =============================================================================#

TEMPLATES = templates_on.h templates_off.h

#-------------------------------

L/a_error.o : a_error.c  config.h gsl_errno.h gsl_message.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/a_stream.o : a_stream.c  config.h gsl_errno.h gsl_message.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/a_results.o : a_results.c  config.h gsl_test.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

#-------------------------------

L/b_init.o : b_init.c  config.h gsl_block.h b_init_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/block.o : block.c  config.h gsl_errno.h gsl_block.h block_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

#-------------------------------

L/c_blas.o : c_blas.c  config.h gsl_math.h gsl_errno.h gsl_cblas.h gsl_blas_types.h gsl_blas.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/c_daxpy.o : c_daxpy.c  gsl_math.h gsl_cblas.h cblas.h c_source_axpy_r.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/c_ddot.o : c_ddot.c    gsl_math.h gsl_cblas.h cblas.h c_source_dot_r.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/c_dnrm2.o : c_dnrm2.c  gsl_math.h gsl_cblas.h cblas.h c_source_nrm2_r.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/c_dscal.o : c_dscal.c  gsl_math.h gsl_cblas.h cblas.h c_source_scal_r.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

#-------------------------------

L/m_directional_minimize.o : m_directional_minimize.c  $(GRADIENT_H)
	$(CC) -o $@ -c  $< $(GSL_INC) 


L/m_conjugate_fr.o : m_conjugate_fr.c  $(GRADIENT_H)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/m_conjugate_pr.o : m_conjugate_pr.c  $(GRADIENT_H)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/m_vector_bfgs.o : m_vector_bfgs.c  $(GRADIENT_H)
	$(CC) -o $@ -c  $< $(GSL_INC) 

#-------------------------------

L/m_convergence.o : m_convergence.c  config.h gsl_multimin.h gsl_blas.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/m_diff.o : m_diff.c  config.h gsl_multimin.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/m_fdfminimizer.o : m_fdfminimizer.c  config.h gsl_errno.h gsl_multimin.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/m_fminimizer.o : m_fminimizer.c config.h gsl_errno.h gsl_multimin.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/m_simplex.o : m_simplex.c  config.h gsl_blas.h gsl_multimin.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/m_steepest_descent.o : m_steepest_descent.c config.h gsl_multimin.h gsl_blas_types.h gsl_blas.h
	$(CC) -o $@ -c  $< $(GSL_INC) 

#-------------------------------
L/v_copy.o : v_copy.c  config.h gsl_vector.h gsl_errno.h v_copy_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/v_init.o : v_init.c  config.h gsl_vector.h v_init_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/v_minmax.o : v_minmax.c  config.h gsl_vector.h v_minmax_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/vector.o : vector.c  config.h gsl_errno.h gsl_vector.h vector_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

#-------------------------------

L/x_getset.o : x_getset.c  config.h gsl_errno.h gsl_matrix.h gsl_vector.h x_getset_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/x_init.o : x_init.c  config.h gsl_errno.h gsl_matrix.h x_init_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

L/x_matrix.o : x_matrix.c  config.h gsl_errno.h gsl_matrix.h x_matrix_source.c $(TEMPLATES)
	$(CC) -o $@ -c  $< $(GSL_INC) 

# =============================================================================#

SLIST = *.c *.h *.sh E/*.sh E/OUT Makefile*  *.cl *.sb TODO* 

ANAME = gslmin
save :
	tar --ignore-failed-read -cvf $(ANAME).tar $(SLIST)
	gzip -6 $(ANAME).tar  
	(cd ../../A/gslmin; mdat) 

# =============================================================================#
