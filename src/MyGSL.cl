;===============================================================================

(load-shared-object "libgslmin.so")
(load-shared-object "libmygsl.so")

(cl:defpackage "TEST-C-CALL" (:use "CL" "SB-ALIEN" "SB-C-CALL"))
(cl:in-package "TEST-C-CALL")

(declaim (sb-ext:muffle-conditions sb-ext:compiler-note))
(declaim (sb-ext:muffle-conditions style-warning))

;(load-shared-object "libgslmin.so")

;//------------------------------------------------------------------------------

(defconstant GSL_SUCCESS    0)
(defconstant GSL_FAILURE   -1)
(defconstant GSL_CONTINUE  -2)

;--------------------------------------------------

(define-alien-type nil (struct gsl-vector

   (size   long)             ; size_t     size;
   (stride long)             ; size_t     stride;
   (data   (* single-float)) ; double    *data;          
   (_block (* long))         ; gsl_block *block; 
   (owner  int)              ; int        owner;

))

;-------------------------------------------------------------------------------
;gsl_vector *gsl_vector_alloc (const size_t n);
;-------------------------------------------------------------------------------
(define-alien-routine gsl-vector-alloc

  (*  (struct gsl-vector)) ; возвращаeмоe значeниe

  (i  long)
)
;-------------------------------------------------------------------------------
;void gsl_vector_set (gsl_vector * v, const size_t i, double x);
;-------------------------------------------------------------------------------
(define-alien-routine gsl-vector-set

  void ; возвращаeмоe значeниe 

  (v  (* (struct gsl-vector)))
  (i  long)
  (x  double)
)
;-------------------------------------------------------------------------------
;double gsl_vector_get (const gsl_vector * v, const size_t i);
;-------------------------------------------------------------------------------
(define-alien-routine gsl-vector-get

  double ; возвращаeмоe значeниe 

  (v  (* (struct gsl-vector)))
  (i  long)
)
;-------------------------------------------------------------------------------

(define-alien-type nil (struct gsl-multimin-fdfminimizer

;  /* multi dimensional part */

  (type  (* (struct gsl_multimin_fdfminimizer_type))) ; const gsl_multimin_fdfminimizer_type *type;
  (fdf   (* (struct gsl_multimin_function_fdf)))      ; gsl_multimin_function_fdf            *fdf;

  (f        double)                   ;  double f;
  (x        (* (struct gsl-vector)))  ;  gsl_vector * x;
  (gradient (* (struct gsl-vector)))  ;  gsl_vector * gradient;
  (dx       (* (struct gsl-vector)))  ;  gsl_vector * dx;

  (state (* t))
))

;-------------------------------------------------------------------------------
;gsl_multimin_fdfminimizer *
;gsl_multimin_fdfminimizer_alloc(const gsl_multimin_fdfminimizer_type *T,
;                                size_t n);
;-------------------------------------------------------------------------------
(define-alien-routine my-gsl-multimin-fdfminimizer-alloc

  (* (struct gsl-multimin-fdfminimizer))              ; возвращаeмоe значeниe 

  ;(type  (* (struct gsl_multimin_fdfminimizer_type))) ; T is already in use 
  (type  c-string)
  (n     long) 
)
;-------------------------------------------------------------------------------
;int 
;gsl_multimin_fdfminimizer_set (gsl_multimin_fdfminimizer * s,
;                               gsl_multimin_function_fdf *fdf,
;                               const gsl_vector * x,
;                               double step_size, double tol);
;-------------------------------------------------------------------------------
(define-alien-routine my-gsl-multimin-fdfminimizer-set

  void ; возвращаeмоe значeниe 

  (s    (* (struct gsl-multimin-fdfminimizer))) 

  (myf   (function double (* (struct gsl-vector)) (* t)))
  (mydf  (function   void (* (struct gsl-vector)) (* t) (* (struct gsl-vector)) ))

  (x      (* (struct gsl-vector)))
  (par    (* (struct gsl-vector)))

  (step_size double)
  (tol       double)
)
;-------------------------------------------------------------------------------
;int
;//gsl_multimin_fdfminimizer_iterate (gsl_multimin_fdfminimizer * s)
;-------------------------------------------------------------------------------
(define-alien-routine gsl-multimin-fdfminimizer-iterate

  ; возвращаeмоe значeниe 
  int 

  (s    (* (struct gsl-multimin-fdfminimizer))) 
)
;-------------------------------------------------------------------------------
;int
;gsl_multimin_test_gradient(const gsl_vector * g, double epsabs);
;-------------------------------------------------------------------------------
(define-alien-routine gsl-multimin-test-gradient

  int  ; возвращаeмоe значeниe 

  (g      (* (struct gsl-vector)))
  (epsabs  double)
)
;-------------------------------------------------------------------------------
;void
;gsl_multimin_fdfminimizer_free (gsl_multimin_fdfminimizer * s)
;-------------------------------------------------------------------------------
(define-alien-routine gsl-multimin-fdfminimizer-free

  void  ; возвращаeмоe значeниe 

  (s    (* (struct gsl-multimin-fdfminimizer))) 
)
;-------------------------------------------------------------------------------
;void
;gsl_vector_free (gsl_vector * v)
;-------------------------------------------------------------------------------
(define-alien-routine gsl-vector-free

  void ; возвращаeмоe значeниe 

  (v  (* (struct gsl-vector)))
)
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------
