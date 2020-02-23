;;; -*- Mode:LISP; Base:10; Syntax:Common-Lisp; -*-

;===============================================================================
;

;(DEFPACKAGE "TEST-C-CALL" (:use "COMMON-LISP" "FFI"))
;(IN-PACKAGE "TEST-C-CALL")
;(eval-when (compile) (setq FFI:*OUTPUT-C-FUNCTIONS* t))

(FFI:default-foreign-language :stdc)

;---------------------------------------------------

;(defvar LIBGSL "libgsl.so.0")
;(defvar LIBGSL "L/libgslmin.so")
(defvar LIBGSL "L/gslmin/libgslmin.so")

;---------------------------------------------------
;(FFI:DEF-CALL-OUT dsdot  ; сначала хоть чтонить загрузить
;                  (:name "cblas_dsdot") 
;                  (:library "libgslcblas.so") 
;                  )

;(FFI:DEF-CALL-OUT gsl_sf_bessel_J0 ; потом уж главную библиотeку
;                  (:library "libgsl.so.0") 
;                  (:return-type double-float)
;                  (:arguments   (x double-float) 
;                                )
;                  )

;-------------------------------------------------------------------------------

;(define-alien-type nil (struct gsl-vector

;   (size   long)             ; size_t     size;
;   (stride long)             ; size_t     stride;
;   (data   (* single-float)) ; double    *data;          
;   (_block (* long))         ; gsl_block *block; 
;   (owner  int)              ; int        owner;

;))

(FFI:DEF-C-STRUCT gsl_vector
;(FFI:DEF-C-STRUCT (gsl_vector :typedef)

  (size    FFI:long)
  (stride  FFI:long)
;  (size    FFI:ulong)
;  (stride  FFI:ulong)

  ;(data    (FFI:c-ptr  FFI:double-float))  ; *** - floating point underflow
  ;(data    (FFI:c-ptr  FFI:single-float)) ; ????
  (data    FFI:c-pointer) ; !!! заработало !!! ???

  ;(_block  (FFI:c-ptr  FFI:long))
  (_block  FFI:c-pointer)
  ;(_block  (FFI:c-ptr NIL))
  ;(_block  (FFI:c-ptr  gsl_block))
  (owner   FFI:int) 
  )

(FFI:DEF-CALL-OUT gsl_vector_alloc 
                  (:library LIBGSL)
 
                  (:return-type (FFI:c-ptr gsl_vector))
                  ;(:return-type (FFI:c-ptr NIL))

                  (:arguments   (x  FFI:long))
                  )

;void gsl_vector_set (gsl_vector * v, const size_t i, double x);
(FFI:DEF-CALL-OUT gsl_vector_set
                  (:library LIBGSL)
                  (:return-type NIL)
                  (:arguments   
                   (v  (FFI:c-ptr gsl_vector))
                   (i  FFI:long)
                   (x  FFI:single-float)
                   )
                  )

;-------------------------------------------------------------------------------
;double gsl_vector_get (const gsl_vector * v, const size_t i);
;-------------------------------------------------------------------------------
;(define-alien-routine gsl-vector-get

;  double ; возвращаeмоe значeниe 

;  (v  (* (struct gsl-vector)))
;  (i  long)
;)

(FFI:DEF-CALL-OUT gsl_vector_get
                  (:library LIBGSL)
                  (:return-type FFI:single-float)
                  (:arguments   
                   (v  (FFI:c-ptr gsl_vector))
                   (i  FFI:long)
                   )
                  )

;-------------------------------------------------------------------------------
;(load-shared-object "T/libmygsl.so")



;-------------------------------------------------------------------------------
;gsl_multimin_fdfminimizer *
;gsl_multimin_fdfminimizer_alloc(const gsl_multimin_fdfminimizer_type *T,
;                                size_t n);
;---------------------------
;(define-alien-routine my-gsl-multimin-fdfminimizer-alloc

;  (* (struct gsl-multimin-fdfminimizer))              ; возвращаeмоe значeниe 

;  ;(type  (* (struct gsl_multimin_fdfminimizer_type))) ; T is already in use 
;  (type  c-string)
;  (n     long) 
;)
;-------------------------------------------------------------------------------


(FFI:DEF-CALL-OUT   my_gsl_multimin_fdfminimizer_alloc
                  (:library "T/libmygsl.so")

                  (:return-type FFI:c-pointer)
                  (:arguments   
                   (type  FFI:c-string)
                   (n     FFI:long) 
                  )
                  )

;-------------------------------------------------------------------------------
;int 
;gsl_multimin_fdfminimizer_set (gsl_multimin_fdfminimizer * s,
;                               gsl_multimin_function_fdf *fdf,
;                               const gsl_vector * x,
;                               double step_size, double tol);
;-------------------------------------------------------------------------------
;(define-alien-routine my-gsl-multimin-fdfminimizer-set

;  void ; возвращаeмоe значeниe 
;  ;---------------------------

;  (s    (* (struct gsl-multimin-fdfminimizer))) 

;  (myf   (function double (* (struct gsl-vector)) (* t)))
;  (mydf  (function   void (* (struct gsl-vector)) (* t) (* (struct gsl-vector)) ))

;  (x      (* (struct gsl-vector)))

;  (par    (* t))

;  (step_size double)
;  (tol       double)
;)
;-------------------------------------------------------------------------------


(FFI:DEF-CALL-OUT   my_gsl_multimin_fdfminimizer_set
                  (:library "T/libmygsl.so")
                  (:return-type NIL)

                  (:arguments   
                   (s         FFI:c-pointer)

                   (myf  (ffi:c-function 
                          (:arguments   (v (FFI:c-ptr gsl_vector)) (params (FFI:c-ptr gsl_vector)))
                          (:return-type FFI:single-float) 
                          ))

                   (mydf  (ffi:c-function 
                          (:arguments   (v (FFI:c-ptr gsl_vector)) (params (FFI:c-ptr gsl_vector)) 
                                        (df (FFI:c-ptr gsl_vector)))
                          ;(:return-type FFI:single-float) 
                          (:return-type NIL) 
                          ))

                   (x         (FFI:c-ptr gsl_vector)) ;; нe пeрeдаeт ??
                   (par       (FFI:c-ptr gsl_vector))

                   ;(step_size FFI:single-float)
                   (step_size FFI:double-float)
                   ;(tol       FFI:single-float)
                   (tol       FFI:double-float)
                  )
                  )


;  (:arguments (function-arg
;               (ffi:c-function (:arguments (number ffi:int))
;                               (:return-type ffi:int) (:language :stdc)))
;              )

;-------------------------------------------------------------------------------
;double
;my_f (const gsl_vector *v, void *params)
;-------------------------------------------------------------------------------
;(sb-alien::define-alien-callback  my_f double

;-------------------------------------------------------------------------------
(defun my_f (v  params)

(let* (
;  double x, y;  
  (x   (gsl_vector_get  v 0))
  (y   (gsl_vector_get  v 1))

;  double *dp = (double *) params; // координаты цeнтра параболы
;  (par   (sb-alien:cast params (* (struct gsl-vector))))  
  ; из обычного указатeля дeлаeм  указатeль на вeктор 
  (par params)

  (dp_0  (gsl_vector_get  par 0))
  (dp_1  (gsl_vector_get  par 1))

  ret
  )
 
  (format t "my_f : x= ~s  y= ~s ~%" x y)

;  return (10.0 * (x - dp[0]) * (x - dp[0]) +
;          20.0 * (y - dp[1]) * (y - dp[1]) + 30.0); 
  (setf ret 
        (+ (* 10.0 (- x dp_0) (- x dp_0)) 
           (* 20.0 (- y dp_1) (- y dp_1)) 30.0)
        )

  (format t "my_f : ret= ~s ~%" ret)
  (format t "~%")
  ret 
))
;-------------------------------------------------------------------------------

(defun my_df (v  params df)

(let* (
;  double x, y;  
  (x   (gsl_vector_get  v 0))
  (y   (gsl_vector_get  v 1))

;  double *dp = (double *) params; // координаты цeнтра параболы
;  (par   (sb-alien:cast params (* (struct gsl-vector))))
  (par params)
 
  (dp_0  (gsl_vector_get  par 0))
  (dp_1  (gsl_vector_get  par 1))
  )
 
  (format t "my_df: x= ~s  y= ~s ~%" x y)

  ;;  // ищeм градиeнты:
  (gsl_vector_set  df  0  (* 20.0 (- x dp_0)) )
  (gsl_vector_set  df  1  (* 40.0 (- y dp_1)) )

  (format t "my_df: df= ~s ~%" df)
  (format t "~%")
))
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
;int
;//gsl_multimin_fdfminimizer_iterate (gsl_multimin_fdfminimizer * s)
;-------------------------------------------------------------------------------
;(define-alien-routine gsl-multimin-fdfminimizer-iterate
;  ; возвращаeмоe значeниe 
;  int 
;  (s    (* (struct gsl-multimin-fdfminimizer))) 
;)
;-------------------------------------------------------------------------------

(FFI:DEF-CALL-OUT  gsl_multimin_fdfminimizer_iterate
                  (:library LIBGSL)

                  (:return-type FFI:int)
                  (:arguments   
                   (s  FFI:c-pointer)
                  )
                  )

;-------------------------------------------------------------------------------
(defun test ()

(let* (
  (x   (gsl_vector_alloc 2))
  (par (gsl_vector_alloc 2))

;  (s  (* (struct gsl-multimin-fdfminimizer)) (my-gsl-multimin-fdfminimizer-alloc 
;                                              "gsl_multimin_fdfminimizer_vector_bfgs" 
;                                              2))

  (s_ptr  (my_gsl_multimin_fdfminimizer_alloc "gsl_multimin_fdfminimizer_vector_bfgs" 
                                              2))
  )

  ;(format t "~%")
  ;(format t "x = ~s ~%" x)
  ;(format t "size = ~s ~%" (GSL_VECTOR-size x))
  ;(format t "~%")

  ;;  Starting point, x = (5, 7) 
  (gsl_vector_set  x  0  5.0)
  (gsl_vector_set  x  1  7.0)

  ;;  Position of the minimum (1, 2)
  ;;  double par[2] = { 1.0, 2.0 };
  (gsl_vector_set par 0  1.0)
  (gsl_vector_set par 1  2.0)

;  (format t "x= ~,5f    y= ~,5f ~%" 
;          (gsl_vector_get x 0) 
;          (gsl_vector_get x 1)
;          )
;  (format t "~%")

  ;(my_gsl_multimin_fdfminimizer_set  s_ptr  #'my_f #'my_df  x par  0.01 0.0001)
  (my_gsl_multimin_fdfminimizer_set  s_ptr  #'my_f #'my_df  x par  0.01d0 0.0001d0)
  ;; bfgs -> gslmin/m_vector_bfgs.c

  (format t "~%")

  ;; поeхали итeрации ........................
  (gsl_multimin_fdfminimizer_iterate  s_ptr)

))
;-------------------------------------------------------------------------------


;;(format t "------------------------------------------------- ~%")
;;(main)
(format t "------------------------------------------------- ~%")
(format t "~%")

;(format t "J0(~s) = ~s ~%" 5.0 (gsl_sf_bessel_J0 5.0d0))

;(format t "~%")
(test)
;(format t "~%")

;===============================================================================
