/* ***BEGIN PROLOGUE DCUHRE */
/* ***DATE WRITTEN   160823   (YYMMDD) */
/* ***REVISION DATE  160823   (YYMMDD) */
/* ***CATEGORY NO.  */
/* ***AUTHOR */
/*            C. C. Kankelborg,, */
/*            Dept. Physics */
/*            Montana State University */
/*                                      */
/*            It was originally written in matlab code, then translated to a c code by SCS */ 
/*            It was modified to find the maximum instead of the minimum for constrained case */   
/*            And Sang-Cheol Seok added a stopping criteion  */
/* ***KEYWORDS Multidimensional Maximization */
/* ***PURPOSE  The routine searches the apporximation of Maximum*/
/*            vector of definite integrals */

/************************************************************
 A simple implementation of the downhill simplex method of
 Nelder & Meade (1963), a.k.a. "amoeba". The function funk
 should be a scalar, real-valued function of a vector argument 
 of the same size as xguess.
 Standard scale values for the reflect function
***************************************************************/


#include "amoeba.h"
#include "utils/utils.h"
#include "config/model.h" // for LINKAGE_DISEQUILIBRIUM

#ifdef STUDYDB
#include "utils/tpl.h"
#include "database/StudyDB.h"
#include "database/databaseSupport.h"
extern struct StudyDB studyDB;
#endif


extern double traitPos;

#define checkpt() fprintf (stderr, "Checkpoint at line %d of file \"%s\"\n",__LINE__,__FILE__)

amoeba_state *held_amoeba_state;


int
amoeba (amoeba_state *sMOD, double xguess[])
{

  int i,j;
  int iterIdx;
  double tempVal;


  /*fprintf (stderr, "The initial gues is eval ");
  for  (j=0; j< sMOD->ndim;j++ ){
    fprintf (stderr, " %f ",xguess[j]);
  }
  fprintf (stderr, " \n");*/
  //checkpt(); 
  /* Step1. Setup the initial guesses based on the current MOD and evaluate at those values*/  
  for  (i=0; i<= sMOD->ndim;i++ ){
    for  (j=0; j< sMOD->ndim;j++ ){
      sMOD->feet[i][j]= xguess[j];
    }
  }

  //checkpt();
  /* Feet includes the initial guess and other n points each of which is only one direction deviated from the initial guess */
  for  (j=0; j< sMOD->ndim;j++ ){
    if((sMOD->xu[j] - xguess[j]) > (xguess[j] - sMOD->xl[j])) {   // xguess[j] is closer to xl[j]
      sMOD->feet[j+1][j] -= (xguess[j]-sMOD->xl[j])/2 ;               
    }else{
      sMOD->feet[j+1][j] += (sMOD->xu[j] - xguess[j])/2 ;
    }
  }
  if (sMOD->verbose >0 ){
    fprintf (stderr, " The sMOD in the beginning of amoeba Just after the initial guess\n");
    MODprint_state(sMOD );
  }

  //checkpt();
  //compute_hlod_2p_qt(xguess, &(tempVal),sMOD->scale);
  //fprintf (stderr, " The sMOD in the beginning of amoeba Just after the initial guess funval =%f\n", tempVal); 
  //checkpt();
  (sMOD->funsub) (xguess, &(tempVal),&(sMOD->scale));
  //fprintf (stderr, " The sMOD in the beginning of amoeba Just after the initial gues funval =%f\n", tempVal);
  //checkpt(); 
  for  (i=0; i<= sMOD->ndim;i++ ){
    
    /*fprintf (stderr, "eval of function at ");
    for  (j=0; j< sMOD->ndim;j++ ){
      fprintf (stderr, " %f ",sMOD->feet[i][j]);
    }
    fprintf (stderr, " \n"); */
    (sMOD->funsub) ( sMOD->feet[i], &(sMOD->funVals[i]),&(sMOD->scale));
    sMOD->total_neval ++;
  }

  sMOD->ztol = sMOD->funVals[0] * 0.001;   /* This is the rel error tolerance for stopping criterion*/
                                     /* sMOD->funVals[0] is the current MOD from DCUHRE */
  sMOD->maxEstMOD = sMOD->funVals[0];

  if (sMOD->verbose >0 ){
    fprintf (stderr, " The sMOD after setting up the feet based on the curr MOD\n");
    MODprint_state(sMOD );
  }

  //checkptsMOD_free(sMOD );();
  /* Step 2 loop for modifying the simplex to improve the max*/
  for (iterIdx =0; sMOD->total_neval  <sMOD->Max_neval; iterIdx++){ 
    pickPQR(sMOD); //sMOD->cur_diff_fail is updated inside pickPQR


  
    if (sMOD->funVals[sMOD->R] > sMOD->maxEstMOD){
      sMOD->maxEstMOD = sMOD->funVals[sMOD->R];
    }

    /* Checking the stopping criterion */
    if ( sMOD->cur_diff_fail >= sMOD->aim_diff_fail){
      break; //return 0;
    }

    reflect_foot(sMOD, sMOD->P , sMOD->reflect);               // Reflect

    if (sMOD->funVals[sMOD->P] > sMOD->funVals[sMOD->R])
      reflect_foot(sMOD, sMOD->P , sMOD->expand);              // Expand;
      
    else if (sMOD->funVals[sMOD->P] < sMOD->funVals[sMOD->Q]){ 
      tempVal = sMOD->funVals[sMOD->P];
      reflect_foot(sMOD, sMOD->P , sMOD->contract);            // 1-dim Contract; 
      if (sMOD->funVals[sMOD->P] > tempVal)
        ndcontract(sMOD, sMOD->R);                             // N-dim contract
    }
    if (sMOD->verbose >0 ){
      fprintf (stderr, " The sMOD after %d iteration. \n", iterIdx);
      MODprint_state(sMOD );
    }
  } 
  
  pickPQR(sMOD);

  if (sMOD->verbose >0 ){
    fprintf (stderr, " The sMOD at the end of amoeba\n");
    MODprint_state(sMOD );
  }

  if (sMOD->funVals[sMOD->R] > sMOD->maxEstMOD){
    sMOD->maxEstMOD = sMOD->funVals[sMOD->R];
  }
  //fprintf (stderr, " The sMOD at the end of amoeba Just after the initial guess funval =%f with %d evals\n", sMOD->maxEstMOD, sMOD->total_neval); 
  //free the structure and allocated memory in initialization
  sMOD_free(sMOD );

  return 0;
}



int 
pickPQR(amoeba_state *sMOD){
  int i;

  /* Identify indices of lowest (P), second lowest (Q), and highest (R) feet. */

  /* Initial guess, certainly wrong */
  sMOD->P=0; 
  sMOD->Q=1; 
  sMOD->R=0; 

  if (sMOD->funVals[sMOD->Q] <  sMOD->funVals[sMOD->P]){ // % Correct the P/Q order for first 2 feet
    sMOD->P=1; 
    sMOD->Q=0;
  }

  for (i=0; i<= sMOD->ndim;i++ ){ //% Loop thru feet, finding P,Q,R
    if (sMOD->funVals[i]  <  sMOD->funVals[sMOD->P]){
      sMOD->Q = sMOD->P; 
      sMOD->P = i;
    }
    else if ((sMOD->funVals[i]  >  sMOD->funVals[sMOD->Q]) && (i != sMOD->P)){
      sMOD->Q = i;
    }
    if (sMOD->funVals[i]  >  sMOD->funVals[sMOD->R])   {
      sMOD->R=i;
    }
  }

  if ( (sMOD->funVals[sMOD->R] - sMOD->maxEstMOD) < sMOD->ztol){
    sMOD->cur_diff_fail ++;
  }
  else{
    sMOD->cur_diff_fail =0;
  }


  return 0;
} //%end pickPQR

int 
reflect_foot(amoeba_state *sMOD, int idx, double scale){
  /* % Reflect the jth foot through the centroid of the other
     % feet of the amoeba. The displacement may be scaled by
     % using scale, whose default value of 1 results in a 
     % reflection that preserves the volume of the amoeba.
     % A scale of 0.5 should never be used, as this would result
     % in a degenerate simplex. Typical scale values:
     %1 ...... reflect through amoeba's opposite face
     %    -0.5 .... stretch the foot outward, doubling amoeba size
     %     0.25 ... shrink inward, halving amoeba size
     % The following variables get updated:
     %     feet(j,:) -- location of the jth foot
     %     f(j) -- value of funk at the jth foot
  **************************************************/
 
  int i,j;
  double cent[sMOD->ndim],disp[sMOD->ndim];
  double noBoundaryPad = 0.001; // This guaranttes the new point is still inside the domain
  //fprintf(stderr," Old point is %f %f %f %f %f\n", sMOD->feet[idx][0],sMOD->feet[idx][1],sMOD->feet[idx][2],sMOD->feet[idx][3],sMOD->feet[idx][4]);   
  memset(cent, 0, sizeof(cent));
  for (j=0; j < sMOD->ndim; j++) {
    for (i=0; i <= sMOD->ndim; i++) {
      cent[j] += sMOD->feet[i][j];
    }
    cent[j] -= sMOD->feet[idx][j];

    cent[j] /= sMOD->ndim;

    disp[j] = 2*(cent[j] - sMOD->feet[idx][j]);

    if (sMOD->feet[idx][j]+scale*disp[j] >= sMOD->xu[j]) {  //check adding 
      sMOD->feet[idx][j] = sMOD->xu[j]-noBoundaryPad ;
    }else if (sMOD->feet[idx][j]+scale*disp[j] <= sMOD->xl[j]) {  //check adding 
      sMOD->feet[idx][j] = sMOD->xl[j]+noBoundaryPad ;
    }else {
      sMOD->feet[idx][j] = sMOD->feet[idx][j] + scale* disp[j];
    }
  } 
  
  //fprintf(stderr,"  idx =%d with scale =%f center is %f %f %f %f %f\n",idx, scale, cent[0],cent[1],cent[2],cent[3],cent[4]);
  //fprintf(stderr,"  disp is %f %f %f %f %f\n", disp[0], disp[1],disp[2],disp[3],disp[4]);
  //fprintf(stderr,"New point is %f %f %f %f %f\n", sMOD->feet[idx][0],sMOD->feet[idx][1],sMOD->feet[idx][2],sMOD->feet[idx][3],sMOD->feet[idx][4]);  
  
  (sMOD->funsub) (sMOD->feet[idx] , &(sMOD->funVals[idx]),&(sMOD->scale));
  sMOD->total_neval ++;
  

  return 0;
}


int 
ndcontract(amoeba_state *sMOD, int idx ){
  /* % Contract all feet, except jth, toward the jth foot.
     % The following variables get updated:
     %     feet -- location of each foot
     %     f -- value of funk at each foot  
  *************************************************************/

  int i,j;

  for (i=0; i <= sMOD->ndim; i++) {
    if (i != idx){
      for (j=0; j < sMOD->ndim; j++) {
        sMOD->feet[i][j] = (sMOD->feet[i][j] + sMOD->feet[idx][j])/2;
      }

      (sMOD->funsub) (sMOD->feet[i] , &(sMOD->funVals[i]),&(sMOD->scale));
      sMOD->total_neval ++;
    }
  }
  
  return 0;
}


int
MODinitialize_state (amoeba_state * sMOD, double *a, double *b, int dim){

  /* This function initializes part of the global variables which can be initialized 
     without any other information.
     More dependent variables are initialized in amoeba_ function after amoeba_
   */
  int i;

  sMOD->ndim = dim;
  sMOD->xl = a;
  sMOD->xu = b;


  sMOD->Max_neval = dim * 100;
  sMOD->maxEstMOD = 0.0;
  sMOD->total_neval = 0;
  sMOD->ifail = 0;

  sMOD->verbose = 0;    /*  0:   no report
                            1:   show how PRQ change over time
                            2:   maximum report  */

  sMOD->reflect =1.0;
  sMOD->expand = -0.5;
  sMOD->contract = 0.25;
 



  sMOD->scale=0;

  /* allocate space for feet and funVals for feet*/
  MALCHOKE(sMOD->feet, sizeof (double *)* (sMOD->ndim+1), double **);
  for (i = 0; i <= sMOD->ndim; i++) {
    MALCHOKE(sMOD->feet[i], sizeof (double) * sMOD->ndim, double *);
  }
  MALCHOKE(sMOD->funVals, sizeof (double) * (sMOD->ndim+1), double *);


  sMOD->cur_diff_fail =0;
  sMOD->aim_diff_fail =dim *2;  // current setting

  return 0;
}



int 
MODprint_state(amoeba_state *sMOD ){

  int i,j;

  fprintf (stderr, "The range is\n");
  for (i=0; i< sMOD->ndim; i++){
    fprintf (stderr, "%d th dimension : [%f , %f] \n",i,sMOD->xl[i],sMOD->xu[i]);
  }
  fprintf (stderr, "The current num of eval = %d with max num of eval = %d\n", sMOD->total_neval, sMOD->Max_neval);

  fprintf (stderr, " Points and funVals with P=%d, Q=%d,  and R=%d\n",sMOD->P,sMOD->Q,sMOD->R);
  for (i=0; i<= sMOD->ndim; i++){
    fprintf (stderr, " %d: ", i);
    for (j=0; j< sMOD->ndim; j++)
      fprintf (stderr, "%f, ",sMOD->feet[i][j]);
    fprintf (stderr, "    %f\n", sMOD->funVals[i]);
  }

  
  fprintf (stderr, "ifail %d\n",sMOD->ifail);
  fprintf (stderr, "verbose %d\n",sMOD->verbose);
  fprintf (stderr, "scale %d\n",sMOD->scale);
  fprintf (stderr, "cur_diff_fail %d\n",sMOD->cur_diff_fail);
  fprintf (stderr, "aim_diff_fail %d\n",sMOD->aim_diff_fail);

  return 0;

}


int 
sMOD_free(amoeba_state *sMOD ){
  int i;

  for (i = 0; i <= sMOD->ndim; i++) {
    free(sMOD->feet[i]);
  }

  free(sMOD->feet);
  free(sMOD->funVals);

  return 0;
}
