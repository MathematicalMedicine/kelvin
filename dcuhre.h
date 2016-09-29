/* ***BEGIN PROLOGUE DCUHRE */
/* ***DATE WRITTEN   900116   (YYMMDD) */
/* ***REVISION DATE  900116   (YYMMDD) */
/* ***CATEGORY NO. H2B1A1 */
/* ***AUTHOR */
/*            Jarle Berntsen, The Computing Centre, */
/*            University of Bergen, Thormohlens gt. 55, */
/*            N-5008 Bergen, Norway */
/*            Phone..  47-5-544055 */
/*            Email..  jarle@eik.ii.uib.no */
/*            Terje O. Espelid, Department of Informatics, */
/*            University of Bergen, Thormohlens gt. 55, */
/*            N-5008 Bergen, Norway */
/*            Phone..  47-5-544180 */
/*            Email..  terje@eik.ii.uib.no */
/*            Alan Genz, Computer Science Department, Washington State */
/*            University, Pullman, WA 99163-2752, USA */
/*            Email..  acg@eecs.wsu.edu */
/* ***KEYWORDS automatic multidimensional integrator, */
/*            n-dimensional hyper-rectangles, */
/*            general purpose, global adaptive */
/* ***PURPOSE  The routine calculates an approximation to a given */
/*            vector of definite integrals */


/***********************************************************************
   Translated in C version by f2c then rewritten to incorporate structures
   by Sang-Cheol Seok
      Center for Quantitative and Computation Biology
      The Research Institute at the Nationwide Children's Hospital
      November 2007.

   Compile and Run

   cc -c test.c
   cc -c dcuhre.c
   cc -o test test.o dcuhre.o -lm
   ./test

***********************************************************************/

#ifndef   __DCUHRE_H__
#define   __DCUHRE_H__


/* header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif




#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef double (*D_fp)(...), (*E_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef double (*D_fp)(), (*E_fp)();
#endif


#define maxdim 30 // 17

enum model_Type {TP_DT=0, TP_QT=1, MP_DT=2, MP_QT=3};

enum ld_Type{LD_ANAL=0, LK_ANAL=1};// Added on 9/28/2016 to accomodate LD case

typedef struct{
  /* all local variables for each sub region kept here */

  int region_id;            /* id for the subregion starting from 0*/
  int parent_id;
  int region_level;         /* level for the subregion  starting from 0*/
  int lchild_id, rchild_id; /* id's for the left and rigt children 
                               0 when the subregion is a leaf*/ 
  double *local_xl, *local_xu; /* local boundary, which is not used currently*/
  double *center, *hwidth;     /*arrays of centers and widths*/


  double local_result;       /* local integraion */
  double local_error;      /* local error estimation */



  /*control variables*/
  int    dir;               /* the directions for further subdivision.*/

  /*Scaling 6/16/2009*/
  int cur_scale;

  int bogusLikelihoods; /* Number of specious likelihoods in this region, if > 0, DO NOT SPLIT */

} sub_region;



typedef struct{
  /* all global variables kept here */
  
  /*Input variables */
  double *xl,*xu;           /* Initial boundary which correspond to a[] and b[] in original code */
  int ndim;                 /* number of variables 1 < NDIM <=  15 */
  int mincls, maxcls;       /* minimum and maximum calls */
  int key;                  /* Key to selected local integration rule. */
  double epsabs, epsrel;    /* Requested absolute and relative errors.*/
  U_fp funsub;              /* function subroutin, which is the integrand */
  int numfun;               /* number of function which is now 1
                               , which should be changed for integrating vector functions*/ 


  /*Output variables */
  double result;            /* final estimation of the integral*/

  double error;             /* final estimation of the error */
  int total_neval;          /* number of total function evaluation. */
  int ifail;                /* indicator of success or reason of fail*/


  /*Control variables*/
  int keyf;                 /* key to selected integration rule. */
  int num;                  /* number of function evaluations over each subregion */
  int wtleng;               /* number of weights in the basic integration rule.  */
  int restar;               /* indicator for the first attempt of integral */  
  int verbose;              /* indicator of running status
                               0 nothing printed
                               1 printing basic subregion information
                               2 printing full sburegio information */
  int maxsub;               /* number of max subregions */
  int minsub;               /* number of min subregions  currently this is not used*/
  int sbrgns;               /* number of current subregeions. Use this for assigning subregion ids*/

  sub_region **sbrg_heap;   /*array of pointers for subregions*/

  /* global variables controling the list of subregions*/
  double greate;             /* Greatest error in the current heap*/
  int next_sbrg;             /* Pointing the subregion with greatest error*/

  /*global variables shared by each subregion, generators, weigts...*/
  double **g;                /* Real array of dimension (NDIM,WTLENG).
                               The fully symmetric sum generators for the rules.
                               G(1,J),...,G(NDIM,J) are the generators for the points
                               associated with the the Jth weights.*/
  double **w;                /* Real array of dimension (5,WTLENG).
                               The weights for the basic and null rules.
                               W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
                               W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.*/
  int *rulpts;              /* the number of points produced by 
                                each generator of the selected rule. */
  double **scales;           /* Scaling factors used to construct new null rules, */
  double **norms;            /* 2**NDIM/(1-norm of the null rule constructed by each of 
                               the scaling factors.) */
  double *errcof;           /* Heuristic error coefficients that are used in the */
                            /* error estimation in BASRUL. */
                            /* It is assumed that the error is computed using: */
                            /* IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3) */
                            /* THEN ERROR = ERRCOF(3)*N1 */
                            /* ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3) */
                            /* ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6)) */
                            /* where N1-N3 are the null rules, EP is the error for */
                            /* the parent */
                            /* subregion and ES is the error for the sibling subregion. */


  /* variables for KELVIN*/
  double *diff_result;         /* diff_error(i)= abs(Br_{i-1} - Br_{i}) Use this for another stopping criterion */
  int nlclass;
  double vol_rate ;    /* Use this to convert results to average function value*/
  enum model_Type mType;
  enum ld_Type ldType;
  
  /* 2-dim arrays for dynamic sampling     3/3/2009 */
  int sampling_mode;       /* turn this on when we apply to Merlin. Default is 0 assigned in initialization fucntion
                           0 : normal dcuhre for Kelvin                             
                           1: sampling mode for one subset
                           2: calculating the BR for only one subset 
                               The main difference is s->maxsub is fixed at one
 */
  double *sample_pts;
  double cur_weight;       /* This holds the weigt for the currently working group of samples */
  int cur_sample;       /* The index of the current sample points*/
  int next_dir;         /* Use this for the next direction to split*/
  /*int sample_dim;         This is always the same as s->ndim   : dimension of dynamic sampling*/
  /*int sample_num;         This is always the same as s->num    : number of sample points*/
  
  /*Scaling 6/16/2009*/
  int scale;
  
  /*Consecutive runs with diff(BR) < error_tol  */
  int cur_diff_suc;
  int aim_diff_suc;  // Aiming number of success, which is set in integrationSupport.c

  /*Consecutive runs with 0<= BR <=0.214  */
  int cur_num_smallBR;
  int aim_num_smallBR;  // Aiming number of success, which is set in integrationSupport.c

} dcuhre_state;





/* Function prototypes */
int initialize_state(dcuhre_state *s,double *, double *, int dim);
int dcuhre_(dcuhre_state *s);
int dadhre_(dcuhre_state *s);
int dinhre_(dcuhre_state *s);
int d132re_(dcuhre_state *s);
int d113re_(dcuhre_state *s);
int d09hre_(dcuhre_state *s);
int d07hre_(dcuhre_state *s);
int dchhre_(dcuhre_state *s);
int drlhre_(dcuhre_state *s, sub_region *cw_sbrg);
int dfshre_(dcuhre_state *s, sub_region *cw_sbrg, double *x, int,double *fulsms, double *funvls);

int print_rule(dcuhre_state *s);
int pow_i(int,int);     /*  integer pow function */
int sbrg_free(sub_region ** sr_control, int);
int print_sbrg(sub_region *cw_sbrg, int dim);





#endif /* __DCUHRE_H__ */
