/* ***BEGIN PROLOGUE DCUHRE */
/* ***DATE WRITTEN   160823   (YYMMDD) */
/* ***REVISION DATE  160823   (YYMMDD) */
/* ***CATEGORY NO.  */
/* ***AUTHOR */
/*            C. C. Kankelborg,, */
/*            Dept. Physics */
/*            Montana State University */
/*                                     */
/*            It was originally written in matlab code, then translated to a c code by SCS */ 
/*            It was modified to find the maximum instead of the minimum for constrained case */   
/*            And Sang-Cheol Seok added a stopping criteion  */
/* ***KEYWORDS Multidimensional Minimization */
/* ***PURPOSE  The routine searches the apporximation of Maximum*/
/*            vector of definite integrals */

/************************************************************
 A simple implementation of the downhill simplex method of
 Nelder & Meade (1963), a.k.a. "amoeba". The function funk
 should be a scalar, real-valued function of a vector argument 
 of the same size as xguess.
 Standard scale values for the reflect function
***************************************************************/


#ifndef   __AMOEBA_H__
#define   __AMOEBA_H__

/* header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dcuhre.h"


typedef struct{

  /*Input variables */
  double *xl,*xu;           /* Initial boundary which correspond to a[] and b[] in original code */
  int ndim;                 /* number of variables 1 < NDIM <=  15 */
  int Max_neval;            /* Max number of evaluation */

  U_fp funsub;              /* function subroutin, which is the integrand */


  /*variables */
  double **feet;  // has ndim+1 points consisting a simplex.
  double *funVals; // has the function values of those points
  
  int P,Q,R; /* Idendify lowest, second lowest, and the highest in feet */

  /* predefined control variables */
  double reflect;
  double expand;
  double contract;


  /*Output variables */
  double maxEstMOD;            /* final estimation of the maximization and used only for the stopping creterion
                                  The estimate is already kept track in integrationSupport.c whenever lk function is called.*/
  int total_neval;          /* number of total function evaluation. */
  int ifail;                /* indicator of success or reason of fail*/

  /*Control variables*/ 
  int verbose;              /* indicator of running status
                               0 nothing printed
                               1 printing basic subregion information
                               2 printing full sburegio information */

  enum model_Type mType;

  /*Scaling 6/16/2009*/
  int scale;


  /*Consecutive runs with the change in maxEst smaller < err_tol  */
  int cur_diff_fail;
  int aim_diff_fail;  // Aiming number of success
  double ztol;       // relative error tolerance for the stopping criterion


} amoeba_state;



/* Function prototype */
int MODinitialize_state(amoeba_state *sMOD, double *, double *, int dim);
int amoeba(amoeba_state *sMOD, double xguess[]);
int pickPQR(amoeba_state *sMOD);
int reflect_foot(amoeba_state *sMOD, int idx, double scale);
int ndcontract(amoeba_state *sMOD, int  idx );

int MODprint_state(amoeba_state *sMOD );
int sMOD_free(amoeba_state *sMOD );

#endif /* __AMOEBA_H__ */
