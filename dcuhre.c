
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
   Translated into C version by f2c then rewritten to incorporate structures
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


#include "dcuhre.h"



#define checkpt() fprintf (stderr, "Checkpoint at line %d of file \"%s\"\n",__LINE__,__FILE__)



int
dcuhre_ (dcuhre_state * s)
{

  int i, num_eval;

  //double fun_val;


  /* Step1.1 Checking inputs  and setting up control variables */
  if (s->verbose > 0) {
    fprintf (stderr, "Checking variables\n");
  }
  dchhre_ (s);


  if (s->ifail != 0) {
    fprintf (stderr, "ifail =%d \n", s->ifail);
    return 1;
  } else {
    if (s->verbose > 0) {
      fprintf (stderr, "Checking variables is done successfully\n");
    }
  }

  /* Step1.2 Setting up global variables controling the list of subregions */
  s->g = (double **) malloc (sizeof (double *) * s->ndim);
  s->w = (double **) malloc (sizeof (double *) * 5);
  s->rulpts = (int *) malloc (sizeof (int) * s->wtleng);
  s->scales = (double **) malloc (sizeof (double *) * 3);
  s->norms = (double **) malloc (sizeof (double *) * 3);
  s->errcof = (double *) malloc (sizeof (double) * 6);
  s->diff_result = (double *) malloc (sizeof (double) * s->maxsub);	/* This is created for another stopping criterion */

  for (i = 0; i < s->ndim; i++) {
    s->g[i] = (double *) malloc (sizeof (double) * s->wtleng);
  }
  for (i = 0; i < 5; i++) {
    s->w[i] = (double *) malloc (sizeof (double) * s->wtleng);
  }
  for (i = 0; i < 3; i++) {
    s->scales[i] = (double *) malloc (sizeof (double) * s->wtleng);
    s->norms[i] = (double *) malloc (sizeof (double) * s->wtleng);
  }



  /* Step1.2 Creating Pointers for subregions */
  s->sbrg_heap = (sub_region **) malloc (sizeof (sub_region *) * s->maxsub);
  for (i = 0; i < s->maxsub; i++) {
    s->sbrg_heap[i] = NULL;
  }


  /* Step2. The adaptive integration routine which calls DTRHRE, DINHRE and DRLHRE */
  if (s->verbose > 0) {
    fprintf (stderr, "Calling dadhre\n");
  }
  s->key = s->keyf;
  dadhre_ (s);

  num_eval = s->total_neval;
  if (s->verbose > 1) {
    print_rule (s);
  }
  //  checkpt();
  //   printf("maxsub is %d  in \n", s->maxsub);
  /* Step 3. free all memory */
  /* free the array of subregions */
  sbrg_free (s->sbrg_heap, s->maxsub);
  //   checkpt();
  free (s->sbrg_heap);
  //    checkpt(); 
  /* free global variables */
  for (i = 0; i < s->ndim; i++) {
    free (s->g[i]);
  }
  for (i = 0; i < 5; i++) {
    free (s->w[i]);
  }
  for (i = 0; i < 3; i++) {
    free (s->scales[i]);
    free (s->norms[i]);
  }

  free (s->g); 
  free (s->w);
  free (s->rulpts);
  free (s->scales);
  free (s->norms);
  free (s->errcof);
  free (s->diff_result);

  return num_eval;
}


int
dadhre_ (dcuhre_state * s)
{
  /*The adaptive integration routine which calls DTRHRE, DINHRE and DRLHRE */
  int i, j, direct;
  int intsgn;			//, pointr;
  sub_region *cw_sbrg, *parent_sbrg;	/* currently working subregion */
  double tmp_result;		/*store the just previous error to calculate the difference in error from the previous one to currrent one */
  double real_result, real_error;

  intsgn = 1;
  i = s->ndim;
  for (j = 0; j < i; ++j) {
    if (s->xu[j] < s->xl[j]) {
      intsgn = -intsgn;
    }
  }


  /* Step2.1 first integraion on the original region
     Call DINHRE to compute the weights and abscissas of */
  /*       the function evaluation points. */
  if (s->verbose > 0) {
    fprintf (stderr, "Initializing generators, weights...\n");
  }
  dinhre_ (s);

  if (s->verbose > 1) {
    fprintf (stderr, "Initial generators, weights, rulpts,..");
    print_rule (s);
  }
  if (s->restar == 1) {
    fprintf (stderr, "This case, restar ==1, is not supported for now\n");
    exit (0);
  }


  /*Step2.2 initialize the first subregion */
  s->sbrg_heap[0] = (sub_region *) malloc (sizeof (sub_region));
  cw_sbrg = s->sbrg_heap[0];

  cw_sbrg->region_id = 0;
  cw_sbrg->region_level = 0;
  cw_sbrg->center = (double *) malloc (sizeof (double) * s->ndim);
  cw_sbrg->hwidth = (double *) malloc (sizeof (double) * s->ndim);
  cw_sbrg->cur_scale = s->scale;
  for (i = 0; i < s->ndim; i++) {
    cw_sbrg->center[i] = (s->xl[i] + s->xu[i]) / 2;
    cw_sbrg->hwidth[i] = fabs ((s->xu[i] - s->xl[i]) / 2);
  }
  if (s->verbose > 1) {
    print_sbrg (cw_sbrg, s->ndim);
  }

  /*Step2.3   Apply DRLHRE over the whole region. */
  if (s->verbose > 0) {
    fprintf (stderr, "Apply DRLHRE over the whole region\n");
  }
  drlhre_ (s, cw_sbrg);
  //fprintf (stderr, "end fo first\n");
  if( cw_sbrg->cur_scale > s->scale){
    s->scale=cw_sbrg->cur_scale;
  }
  s->total_neval += s->num;
  s->result += cw_sbrg->local_result;
  s->error += cw_sbrg->local_error;
  tmp_result = s->result /s->vol_rate;


  if (s->verbose > 0) {
    fprintf (stderr, "After %d sub regions result=%20.18f error=%20.18f\n",
	     s->sbrgns, s->result, s->error);
  }


  s->greate = cw_sbrg->local_error;
  s->sbrgns++;
  s->next_sbrg = 0;
  /*Some action which DRTHRE_ does should be here */


  /* ***End initialisation. */
  if(s->sampling_mode){  /* Only one subregion for sampling */
    if(s->sampling_mode==2)
      s->next_dir = cw_sbrg->dir;
    return 0;
  }

  /*Step3 Loop for main integration */
  while (s->sbrgns < s->maxsub) {

    /*This is for setting absolute error */
    //if(s->mType == MT_DT){
    s->epsabs = 0.01;		/*This is error tolerance */

    real_result = s->result / s->vol_rate; 
    real_error = s->error / s->vol_rate;

    if (s->verbose > 0) 
      fprintf(stderr, "s result %f  real result %f\n ", s->result, real_result);

    s->epsabs *=
      (-5.77 + 54.0 * real_result + real_result * real_result) * (-5.77 +
								  54.0 *
								  real_result
								  +
								  real_result
								  *
								  real_result);
    s->epsabs /= (-11.54 * real_result + 54.0 * real_result * real_result);

    if(s->epsabs <0)
      s->epsabs =0.0;


    //s->epsrel = s->epsabs / (s->result);
    if (s->sbrgns == 1) {
      s->diff_result[0] = s->epsabs * 2;	/*Dummy number for the main while loop */
    }
    if (s->verbose > 0) 
      fprintf (stderr, "Setting absolute error %12.8f  diff =%f  real error = %f real result=%f \n",
	     s->epsabs, s->diff_result[s->sbrgns - 1], real_error, real_result);


    //}

    //if ( s->error > s->epsabs) {    //this is for stopping criterion 
    //if ((s->error > s->epsrel *fabs(s->result)) && (s->error > s->epsabs)) {    
    //if ((s->diff_result[s->sbrgns - 1] > s->epsabs) && (s->error > s->epsabs)) { // before 5/18/2008
    //if ((s->result <0)||(s->diff_result[s->sbrgns - 1] > s->epsabs) || real_result <0.25) {
    //if ((s->result <0)||((s->diff_result[s->sbrgns - 1] > s->epsabs) && (real_error > s->epsabs))) {  // before 11/25/2008
    if (((s->result <0)|| (s->diff_result[s->sbrgns - 1] >= s->epsabs)) && ( (real_result <0.9)|| (real_error > s->epsabs) )){// adding s->result <0) on 3/23/2009
      /*   If we are allowed to divide further, */
      /*   prepare to apply basic rule over each half of the */
      /*   NDIV subregions with greatest errors. */
      /*   If MAXSUB is great enough, NDIV = MDIV */
      if (s->verbose > 0) {
	fprintf (stderr,
		 "Select %d subregion with %f error and direction = %d\n",
		 s->next_sbrg, s->greate, s->sbrg_heap[s->next_sbrg]->dir);
      }
      /*Step 3.1   Adjust RESULT and ABSERR. and prepare for adding children */
      parent_sbrg = s->sbrg_heap[s->next_sbrg];	/* parent_sbrg is the sbrg to be split */

      s->result -= parent_sbrg->local_result;
      s->error -= parent_sbrg->local_error;

      /*Step 3.2   Compute first half region (left child) */
      s->sbrg_heap[s->sbrgns] = (sub_region *) malloc (sizeof (sub_region));
      cw_sbrg = s->sbrg_heap[s->sbrgns];

      cw_sbrg->parent_id = parent_sbrg->region_id;
      cw_sbrg->region_id = s->sbrgns;
      parent_sbrg->lchild_id = cw_sbrg->region_id;
      cw_sbrg->region_level = parent_sbrg->region_level + 1;
      cw_sbrg->center = (double *) malloc (sizeof (double) * s->ndim);
      cw_sbrg->hwidth = (double *) malloc (sizeof (double) * s->ndim);
      cw_sbrg->cur_scale = s->scale;
      for (i = 0; i < s->ndim; i++) {
	cw_sbrg->center[i] = parent_sbrg->center[i];
	cw_sbrg->hwidth[i] = parent_sbrg->hwidth[i];
      }
      direct = parent_sbrg->dir;
      cw_sbrg->hwidth[direct] = parent_sbrg->hwidth[direct] / 2;
      cw_sbrg->center[direct] += cw_sbrg->hwidth[direct];
      cw_sbrg->lchild_id = 0;

      if (s->verbose > 1) {
	print_sbrg (cw_sbrg, s->ndim);
      }
      drlhre_ (s, cw_sbrg);
      if( cw_sbrg->cur_scale > s->scale){
        s->scale=cw_sbrg->cur_scale;
      }
      s->total_neval += s->num;

      s->result += cw_sbrg->local_result;
      s->error += cw_sbrg->local_error;
      if (s->verbose > 1) {
	fprintf (stderr, "After %d sub regions result=%10.8f error=%10.8f\n",
		 s->sbrgns, s->result, s->error);
      }
      s->sbrgns++;

      /*Step 3.3 Compute second half region.(right child) */
      s->sbrg_heap[s->sbrgns] = (sub_region *) malloc (sizeof (sub_region));
      cw_sbrg = s->sbrg_heap[s->sbrgns];

      cw_sbrg->parent_id = parent_sbrg->region_id;
      cw_sbrg->region_id = s->sbrgns;
      parent_sbrg->rchild_id = cw_sbrg->region_id;
      cw_sbrg->region_level = parent_sbrg->region_level + 1;
      cw_sbrg->center = (double *) malloc (sizeof (double) * s->ndim);
      cw_sbrg->hwidth = (double *) malloc (sizeof (double) * s->ndim);
      cw_sbrg->cur_scale = s->scale;
      for (i = 0; i < s->ndim; i++) {
	cw_sbrg->center[i] = parent_sbrg->center[i];
	cw_sbrg->hwidth[i] = parent_sbrg->hwidth[i];
      }

      cw_sbrg->hwidth[direct] = parent_sbrg->hwidth[direct] / 2;
      cw_sbrg->center[direct] -= cw_sbrg->hwidth[direct];
      cw_sbrg->lchild_id = 0;

      if (s->verbose > 1) {
	print_sbrg (cw_sbrg, s->ndim);
      }
      drlhre_ (s, cw_sbrg);
      if( cw_sbrg->cur_scale > s->scale){
        s->scale=cw_sbrg->cur_scale;
      }
      s->total_neval += s->num;

      s->result += cw_sbrg->local_result;
      s->error += cw_sbrg->local_error;
      real_result = s->result / s->vol_rate;
      real_error = s->error / s->vol_rate;

      s->diff_result[s->sbrgns] = fabs (real_result  - tmp_result);
      tmp_result = real_result;
      if (s->verbose > 0) {
	fprintf (stderr,
		 "After %d sub regions result=%10.8f error=%10.8f and difference=%10.8f\n",
		 s->sbrgns, real_result, real_error, s->diff_result[s->sbrgns]);
      }

      s->sbrgns++;


      /*find the next subregion with greatest error among leaves to split */
      s->greate = 0.0;
      for (i = 0; i < s->sbrgns; i++) {
	cw_sbrg = s->sbrg_heap[i];
	if (cw_sbrg->lchild_id == 0) {	/*That is this sbrg is a leaf */
	  if (cw_sbrg->local_error > s->greate) {
	    s->next_sbrg = cw_sbrg->region_id;
	    s->greate = cw_sbrg->local_error;
	  }
	}
      }
      s->ifail = 1;
    } else {
      s->ifail = 0;
      break;
    }
    
    if((s->sbrgns  == s->maxsub) && (s->result / s->vol_rate <0)){   
                  /* in case of a negative result, double the current maxsub*/
      s->maxsub = 2* s->maxsub+1;
      s->diff_result=realloc(s->diff_result, sizeof (double) *s->maxsub);
      s->sbrg_heap=realloc(s->sbrg_heap, sizeof (sub_region *) * s->maxsub);
      for (i = s->maxsub/2; i < s->maxsub; i++) {
        s->sbrg_heap[i] = NULL;
      }
      //      printf("maxsub is doubled to %d \n",s->maxsub );
    } 
     
  }
  //  checkpt();
  //  printf("maxsub is %d in dadddd\n", s->maxsub);

  return 0;

}


int
drlhre_ (dcuhre_state * s, sub_region * cw_sbrg)
{

  int i, k, divaxn;
  double rgnvol, difsum, difmax, frthdf, ratio, search;
  double d__1, d__2, d__3;
  double *x, *null;
  int g_work_col;

  rgnvol = 1.0;
  divaxn = 0;

  x = (double *) malloc (sizeof (double) * s->ndim);
  null = (double *) malloc (sizeof (double) * 8);
  for (i = 0; i < s->ndim; i++) {
    rgnvol *= cw_sbrg->hwidth[i];
    x[i] = cw_sbrg->center[i];
    if (cw_sbrg->hwidth[i] > cw_sbrg->hwidth[divaxn]) {
      divaxn = i;
    }
  }

  if(s->sampling_mode>0)
    (s->funsub) ( x, &(cw_sbrg->local_error),s);
  else{
    (s->funsub) ( x, &(cw_sbrg->local_error), &(cw_sbrg->cur_scale));
  }
  if(s->sampling_mode==1)
    s->sample_pts[s->cur_sample*(s->ndim+1)] = s->w[0][0];

  //fprintf (stderr, "1 local result =%10.8f  and local error =%10.8f\n", cw_sbrg->local_result,cw_sbrg->local_error);
  //fprintf(stderr," sw00 is %G   %p %p\n", s->w[0][0], s->w, s->w[0]);

  if(isnan(s->w[0][0])){
    fprintf(stderr,"sw00 is nan \n");
    exit(1);
  }

  cw_sbrg->local_result = s->w[0][0] * cw_sbrg->local_error;	//rgnerr[j];

  for (k = 0; k < 4; ++k) {
    null[k] = s->w[k + 1][0] * cw_sbrg->local_error;	//  rgnerr[j];
  }

  difmax = 0.;
  d__1 = s->g[0][2] / s->g[0][1];
  ratio = d__1 * d__1;
  for (i = 0; i < s->ndim; i++) {
    x[i] = cw_sbrg->center[i] - cw_sbrg->hwidth[i] * s->g[0][1];

    if(s->sampling_mode>0)
      (s->funsub) ( x, &null[4],s);
    else{
      (s->funsub) ( x, &null[4], &(cw_sbrg->cur_scale));
    }
    if(s->sampling_mode==1)
      s->sample_pts[s->cur_sample*(s->ndim+1)] = s->w[0][1];

    x[i] = cw_sbrg->center[i] + cw_sbrg->hwidth[i] * s->g[0][1];

    if(s->sampling_mode>0)
      (s->funsub) ( x, &null[5],s);
    else{
      (s->funsub) ( x, &null[5], &(cw_sbrg->cur_scale));
    }
    if(s->sampling_mode==1)
      s->sample_pts[s->cur_sample*(s->ndim+1)] =  s->w[0][1];

    x[i] = cw_sbrg->center[i] - cw_sbrg->hwidth[i] * s->g[0][2];
    if(s->sampling_mode>0)
      (s->funsub) ( x, &null[6],s);
    else{
      (s->funsub) ( x, &null[6], &(cw_sbrg->cur_scale));
    }
    if(s->sampling_mode==1)
      s->sample_pts[s->cur_sample*(s->ndim+1)] = s->w[0][2];

    x[i] = cw_sbrg->center[i] + cw_sbrg->hwidth[i] * s->g[0][2];
    if(s->sampling_mode>0)
      (s->funsub) (x, &null[7],s);
    else{
      (s->funsub) ( x, &null[7], &(cw_sbrg->cur_scale));
    }
    if(s->sampling_mode==1)
      s->sample_pts[s->cur_sample*(s->ndim+1)] = s->w[0][2];

    x[i] = cw_sbrg->center[i];
    difsum = 0.;

    frthdf = (1 - ratio) * 2 * cw_sbrg->local_error - (null[6] +
						       null[7]) +
      ratio * (null[4] + null[5]);

    //fprintf (stderr, "printing null 4 %f %f %f %f %f\n", null[4],null[5],null[6],null[7],frthdf);
    /* Ignore differences below roundoff */
    if (cw_sbrg->local_error + frthdf / 4 != cw_sbrg->local_error) {
      difsum += fabs (frthdf);
    }
    for (k = 0; k < 4; k++) {
      null[k] += s->w[k + 1][1] * (null[4] + null[5])
	+ s->w[k + 1][2] * (null[6] + null[7]);
    }

    cw_sbrg->local_result +=
      s->w[0][1] * (null[4] + null[5]) + s->w[0][2] * (null[6] + null[7]);

    if (difsum > difmax) {
      difmax = difsum;
      divaxn = i;
    }
  }
  cw_sbrg->dir = divaxn;

  //fprintf (stderr, "null is %10.8f  %10.8f  %10.8f \n", null[0], null[1], null[2]);
  /* Finish computing the rule values. */
  for (i = 3; i < s->wtleng; i++) {
    //dfshre_(s->ndim, cw_sbrg->center, cw_sbrg->hwidth, x, s->g[0][i], 
    //                s->numfun, (S_fp)s->funsub, &(cw_sbrg->local_error), &null[4]);
    g_work_col = i;
    //fprintf (stderr, "i=%d calling dfshre\n",i);
    //fprintf (stderr, "g is %f %f %f %f %f\n",s->g[0][g_work_col],s->g[1][g_work_col],s->g[2][g_work_col],s->g[3][g_work_col],s->g[4][g_work_col]);

    if(s->sampling_mode==1){
      s->cur_weight = s->w[0][i];
    }


    dfshre_ (s, cw_sbrg, x, g_work_col, &(cw_sbrg->local_error), &null[4]);
    cw_sbrg->local_result += s->w[0][i] * cw_sbrg->local_error;

    for (k = 0; k < 4; k++) {
      null[k] += s->w[k + 1][i] * cw_sbrg->local_error;
    }

  }

  if (s->verbose > 0) {

    fprintf (stderr, "local result =%10.8f  and local error =%10.8f\n", cw_sbrg->local_result,cw_sbrg->local_error);
  }
  if(isnan(cw_sbrg->local_result)){
    fprintf(stderr,"local result is nan \n");
    exit(1);
  }

  /*    Compute errors. */
  /*    We search for the null rule, in the linear space spanned by two */
  /*    successive null rules in our sequence, which gives the greatest */
  /*    error estimate among all normalized (1-norm) null rules in this */
  /*    space. */
  //fprintf (stderr, "null is %10.8f  %10.8f  %10.8f \n", null[0], null[1], null[2]);
  for (i = 0; i < 3; i++) {
    search = 0.;

    for (k = 0; k < s->wtleng; k++) {
      /* Computing MAX */
      d__2 = search, d__3 = (d__1 = null[i + 1] + s->scales[i][k] * null[i],
			     fabs (d__1)) * s->norms[i][k];
      search = max (d__2, d__3);
    }
    null[i] = search;
  }
  //fprintf (stderr, "null is %10.8f  %10.8f  %10.8f \n", null[0], null[1], null[2]);
  if (s->errcof[0] * null[0] <= null[1] && s->errcof[1] * null[1] <= null[2]) {
    cw_sbrg->local_error = s->errcof[2] * null[0];
  } else {
    /* Computing MAX */
    d__1 = null[0], d__2 = null[1], d__1 = max (d__1, d__2), d__2 = null[2];
    cw_sbrg->local_error = s->errcof[3] * max (d__1, d__2);
  }
  cw_sbrg->local_error *= rgnvol;
  cw_sbrg->local_result *= rgnvol;

  //fprintf (stderr, "local result =%10.8f  and local error =%10.8f\n", cw_sbrg->local_result,cw_sbrg->local_error);

  free (x);
  free (null);

  return 0;
}


int dfshre_ (dcuhre_state * s, sub_region * cw_sbrg, double *x, int g_work_col,double *fulsms, double *funvls)
{

  /* Local variables */
  int i, l, i__1, i__2;
  double gi, gl;
  int ixchng = 0, lxchng = 0;

  *fulsms = 0.;

  // fprintf(stderr," sw00 is %G   %p %p\n", s->w[0][0], s->w, s->w[0]);
  /* Compute centrally symmetric sum for permutation of G */
L20:

  for (i = 0; i < s->ndim; i++) {
    x[i] = cw_sbrg->center[i] + s->g[i][g_work_col] * cw_sbrg->hwidth[i];
  }

L40:
  if(s->sampling_mode>0)  
    (s->funsub) ( x, funvls,s);
  else
    (s->funsub) ( x, funvls, &(cw_sbrg->cur_scale));

  if(s->sampling_mode==1)
    s->sample_pts[s->cur_sample*(s->ndim+1)] = s->cur_weight;

  *fulsms += *funvls;

  for (i = 0; i < s->ndim; i++) {
    if (fabs (s->g[i][g_work_col]) > 0.000000000001)
      s->g[i][g_work_col] = -s->g[i][g_work_col];
    x[i] = cw_sbrg->center[i] + s->g[i][g_work_col] * cw_sbrg->hwidth[i];
    if (s->g[i][g_work_col] < 0.) {
      goto L40;
    }
  }
  //fprintf (stderr, "go %f %f %f %f %f\n",s->g[0][g_work_col],s->g[1][g_work_col],s->g[2][g_work_col],s->g[3][g_work_col],s->g[4][g_work_col]);
  /*   Find next distinct permutation of G and loop back for next sum. */
  /*   Permutations are generated in reverse lexicographic order. */
  for (i = 2; i <= s->ndim; i++) {
    //fprintf (stderr, "i=%d g= %f   g=%f\n",i,s->g[i - 2][g_work_col] , s->g[i-1][g_work_col]);
    if (s->g[i - 2][g_work_col] > s->g[i - 1][g_work_col]) {
      gi = s->g[i - 1][g_work_col];
      ixchng = i - 1;
      i__2 = (i - 1) / 2;
      for (l = 1; l <= i__2; l++) {
	gl = s->g[l - 1][g_work_col];
	s->g[l - 1][g_work_col] = s->g[i - l - 1][g_work_col];
	s->g[i - l - 1][g_work_col] = gl;
	if (gl <= gi) {
	  --ixchng;
	}
	if (s->g[l - 1][g_work_col] > gi) {
	  //fprintf (stderr, "lxchng =%d is changing to l=%d\n",lxchng,l);
	  lxchng = l;
	}
	// fprintf (stderr, "     i=%d  l=%d  IXCHNG =%d, lxchng =%d\n", i, l, ixchng, lxchng);
      }
      if (s->g[ixchng - 1][g_work_col] <= gi) {
	ixchng = lxchng;
      }
      s->g[i - 1][g_work_col] = s->g[ixchng - 1][g_work_col];
      s->g[ixchng - 1][g_work_col] = gi;
      //fprintf (stderr, "gf %f %f %f %f %f IXCHNG =%d, lxchng =%d\n",s->g[0][g_work_col],s->g[1][g_work_col],s->g[2][g_work_col],s->g[3][g_work_col],s->g[4][g_work_col], ixchng, lxchng);
      goto L20;
    }

/* L80: */
  }

/*     Restore original order to generators */
  //checkpt();
  i__1 = s->ndim / 2;
  for (i = 1; i <= i__1; i++) {
    gi = s->g[i - 1][g_work_col];
    s->g[i - 1][g_work_col] = s->g[s->ndim - i][g_work_col];
    s->g[s->ndim - i][g_work_col] = gi;
  }

  return 0;
}				/* End of dfshre_ */




int
dinhre_ (dcuhre_state * s)
{

  int i, j, k;
  double we[14];		/*temporary array */

  if (s->key == 1) {
    d132re_ (s);
  } else if (s->key == 2) {
    d113re_ (s);
  } else if (s->key == 3) {
    d09hre_ (s);
  } else if (s->key == 4) {
    d07hre_ (s);
  }

  /*   Compute SCALES and NORMS. */
  for (k = 0; k < 3; k++) {
    for (i = 0; i < s->wtleng; i++) {
      if (s->w[k + 1][i] != 0.) {
	s->scales[k][i] = -s->w[k + 2][i] / s->w[k + 1][i];
      } else {
	s->scales[k][i] = 100.;
      }

      for (j = 0; j < s->wtleng; j++) {
	we[j] = s->w[k + 2][j] + s->scales[k][i] * s->w[k + 1][j];
      }
      s->norms[k][i] = 0.;
      for (j = 0; j < s->wtleng; j++) {
	s->norms[k][i] += s->rulpts[j] * fabs (we[j]);
      }
      s->norms[k][i] = pow_i (2, s->ndim) / s->norms[k][i];
    }
  }

  return 0;

}

int
d132re_ (dcuhre_state * s)
{
  /* Initialized data */

  double dim2g[16] = { .2517129343453109, .7013933644534266,
    .9590960631619962, .9956010478552127, .5, .1594544658297559,
    .3808991135940188, .6582769255267192, .8761473165029315,
    .998243184053198, .9790222658168462, .6492284325645389,
    .8727421201131239, .3582614645881228, .5666666666666666,
    .2077777777777778
  };
  double dim2w[70] /* was [14][5] */  = { .0337969236013446,
    .09508589607597761, .1176006468056962, .0265777458632695,
    .0170144177020064, 0., .0162659309863741, .1344892658526199,
    .1328032165460149, .0563747476999187, .0039082790813105,
    .0301279877743215, .1030873234689166, .0625, .3213775489050763,
    -.1767341636743844, .07347600537466072, -.03638022004364754,
    .02125297922098712, .1460984204026913, .01747613286152099,
    .1444954045641582, 1.307687976001325e-4, 5.380992313941161e-4,
    1.042259576889814e-4, -.001401152865045733, .008041788181514763,
    -.1420416552759383, .3372900883288987, -.1644903060344491,
    .07707849911634622, -.0380447835850631, .02223559940380806,
    .1480693879765931, 4.467143702185814e-6, .150894476707413,
    3.647200107516215e-5, 5.77719899901388e-4, 1.041757313688177e-4,
    -.001452822267047819, .008338339968783705, -.147279632923196,
    -.8264123822525677, .306583861409436, .002389292538329435,
    -.1343024157997222, .088333668405339, 0., 9.786283074168292e-4,
    -.1319227889147519, .00799001220015063, .003391747079760626,
    .002294915718283264, -.01358584986119197, .04025866859057809,
    .003760268580063992, .6539094339575232, -.2041614154424632,
    -.174698151579499, .03937939671417803, .006974520545933992, 0.,
    .006667702171778258, .05512960621544304, .05443846381278607,
    .02310903863953934, .01506937747477189, -.0605702164890189,
    .04225737654686337, .02561989142123099
  };

  int i, j;

/*   Assign values to W. */
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 14; j++) {
      s->w[i][j] = dim2w[i * 14 + j];
    }
  }

/*   Assign values to G. */
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 14; ++j) {
      s->g[i][j] = 0.;
    }
  }
  s->g[0][1] = dim2g[0];
  s->g[0][2] = dim2g[1];
  s->g[0][3] = dim2g[2];
  s->g[0][4] = dim2g[3];
  s->g[0][5] = dim2g[4];
  s->g[0][6] = dim2g[5];
  s->g[1][6] = s->g[0][6];
  s->g[0][7] = dim2g[6];
  s->g[1][7] = s->g[0][7];
  s->g[0][8] = dim2g[7];
  s->g[1][8] = s->g[0][8];
  s->g[0][9] = dim2g[8];
  s->g[1][9] = s->g[0][9];
  s->g[0][10] = dim2g[9];
  s->g[1][10] = s->g[0][10];
  s->g[0][11] = dim2g[10];
  s->g[1][11] = dim2g[11];
  s->g[0][12] = dim2g[12];
  s->g[1][12] = dim2g[13];
  s->g[0][13] = dim2g[14];
  s->g[1][13] = dim2g[15];

/*   Assign values to RULPTS. */

  s->rulpts[0] = 1.;
  for (i = 1; i < 11; i++) {
    s->rulpts[i] = 4.;
  }
  s->rulpts[11] = 8.;
  s->rulpts[12] = 8.;
  s->rulpts[13] = 8.;

/*   Assign values to ERRCOF. */

  s->errcof[0] = 10.;
  s->errcof[1] = 10.;
  s->errcof[2] = 1.f;
  s->errcof[3] = 5.f;
  s->errcof[4] = .5f;
  s->errcof[5] = .25f;


  return 0;
}				/* End of d132re_ */


int
d113re_ (dcuhre_state * s)
{

  /* Initialized data */

  double dim3g[14] = { .19, .5, .75, .8, .9949999999999999,
    .99873449983514, .7793703685672423, .9999698993088767,
    .7902637224771788, .4403396687650737, .4378478459006862,
    .9549373822794593, .9661093133630748, .4577105877763134
  };
  double dim3w[65] /* was [13][5] */  = { .007923078151105734,
    .0679717739278808, .001086986538805825, .1838633662212829,
    .03362119777829031, .01013751123334062, .001687648683985235,
    .1346468564512807, .001750145884600386, .07752336383837454,
    .2461864902770251, .06797944868483039, .01419962823300713,
    1.715006248224684, -.3755893815889209, .1488632145140549,
    -.2497046640620823, .1792501419135204, .00344612675897389,
    -.005140483185555825, .006536017839876425, -6.5134549392297e-4,
    -.006304672433547204, .01266959399788263, -.005454241018647931,
    .004826995274768427, 1.936014978949526, -.3673449403754268,
    .02929778657898176, -.1151883520260315, .05086658220872218,
    .04453911087786469, -.022878282571259, .02908926216345833,
    -.002898884350669207, -.02805963413307495, .05638741361145884,
    -.02427469611942451, .02148307034182882, .517082819560576,
    .01445269144914044, -.3601489663995932, .3628307003418485,
    .007148802650872729, -.09222852896022966, .01719339732471725,
    -.102141653746035, -.007504397861080493, .01648362537726711,
    .05234610158469334, .01445432331613066, .003019236275367777,
    2.05440450381852, .0137775998849012, -.576806291790441,
    .03726835047700328, .006814878939777219, .05723169733851849,
    -.04493018743811285, .02729236573866348, 3.54747395055699e-4,
    .01571366799739551, .04990099219278567, .0137791555266677,
    .002878206423099872
  };

  int i, j;



/* ***FIRST EXECUTABLE STATEMENT D113RE */

/*   Assign values to W. */
  for (i = 0; i < 5; i++) {
    for (j = 0; j < 13; j++) {
      s->w[i][j] = dim3w[i * 13 + j];
    }
  }

/*   Assign values to G. */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 13; j++) {
      s->g[i][j] = 0.;
    }
  }
  s->g[0][1] = dim3g[0];
  s->g[0][2] = dim3g[1];
  s->g[0][3] = dim3g[2];
  s->g[0][4] = dim3g[3];
  s->g[0][5] = dim3g[4];
  s->g[0][6] = dim3g[5];
  s->g[1][6] = s->g[0][6];
  s->g[0][7] = dim3g[6];
  s->g[1][7] = s->g[0][7];
  s->g[0][8] = dim3g[7];
  s->g[1][8] = s->g[0][8];
  s->g[2][8] = s->g[0][8];
  s->g[0][9] = dim3g[8];
  s->g[1][9] = s->g[0][9];
  s->g[2][9] = s->g[0][9];
  s->g[0][10] = dim3g[9];
  s->g[1][10] = s->g[0][10];
  s->g[2][10] = s->g[0][10];
  s->g[0][11] = dim3g[11];
  s->g[1][11] = dim3g[10];
  s->g[2][11] = s->g[1][11];
  s->g[0][12] = dim3g[12];
  s->g[1][12] = s->g[0][12];
  s->g[2][12] = dim3g[13];

/*   Assign values to RULPTS. */

  s->rulpts[0] = 1.;
  s->rulpts[1] = 6.;
  s->rulpts[2] = 6.;
  s->rulpts[3] = 6.;
  s->rulpts[4] = 6.;
  s->rulpts[5] = 6.;
  s->rulpts[6] = 12.;
  s->rulpts[7] = 12.;
  s->rulpts[8] = 8.;
  s->rulpts[9] = 8.;
  s->rulpts[10] = 8.;
  s->rulpts[11] = 24.;
  s->rulpts[12] = 24.;

/*   Assign values to ERRCOF. */

  s->errcof[0] = 4.;
  s->errcof[1] = 4.f;
  s->errcof[2] = .5f;
  s->errcof[3] = 3.f;
  s->errcof[4] = .5f;
  s->errcof[5] = .25f;


  return 0;
}				/* End of d113re_ */



int
d09hre_ (dcuhre_state * s)
{

  double d;

  /* Local variables */
  int i, j;
  double lam0, lam1, lam2, lam3, lamp, ratio, twondm;

  for (j = 0; j < s->wtleng; j++) {
    for (i = 0; i < s->ndim; i++) {
      s->g[i][j] = 0.0;
    }
    for (i = 0; i < 5; i++) {
      s->w[i][j] = 0.0;
    }
    s->rulpts[j] = 2 * s->ndim;
  }

  twondm = (double) pow_i (2, s->ndim);
  s->rulpts[s->wtleng - 1] = twondm;
  if (s->ndim > 2)
    s->rulpts[7] = (double) (s->ndim * 4 * (s->ndim - 1) * (s->ndim - 2) / 3);

  s->rulpts[6] = 4. * s->ndim * (s->ndim - 1);
  s->rulpts[5] = 2. * s->ndim * (s->ndim - 1);
  s->rulpts[0] = 1.0;

  /*     Compute squared generator parameters */
  lam0 = .4707f;
  lam1 = 4 / (15 - 5 / lam0);
  ratio = (1 - lam1 / lam0) / 27;
  lam2 =
    (5 - lam1 * 7 - ratio * 35) / (7 - lam1 * 35 / 3 - ratio * 35 / lam0);
  ratio = ratio * (1 - lam2 / lam0) / 3;
  lam3 =
    (7 - (lam2 + lam1) * 9 + lam2 * 63 * lam1 / 5 - ratio * 63) / (9 -
								   (lam2 +
								    lam1) *
								   63 / 5 +
								   lam2 * 21 *
								   lam1 -
								   ratio *
								   63 / lam0);
  lamp = .0625f;

  //fprintf (stderr, "lam0 = %10.7f  %10.7f %10.7f %10.7f %10.7f \n", lam0,lam1,lam2,lam3,lamp);

/*     Compute degree 9 rule weights */
  d = 3 * lam0;
  d *= d * d * d;		/* (3* lam0)^4 */
  s->w[0][s->wtleng - 1] = 1 / d / twondm;

  if (s->ndim > 2) {
    d = 6 * lam1;
    d *= d * d;			/* (6* lamd1)^3 */
    s->w[0][7] = (1 - 1 / (lam0 * 3)) / d;
  }

  s->w[0][6] =
    (1 - (lam0 + lam1) * 7 / 5 +
     lam0 * 7 * lam1 / 3) / (lam1 * 84 * lam2 * (lam2 - lam0) * (lam2 -
								 lam1));
  s->w[0][5] =
    (1 - (lam0 + lam2) * 7 / 5 +
     lam0 * 7 * lam2 / 3) / (lam1 * 84 * lam1 * (lam1 - lam0) * (lam1 -
								 lam2)) -
    s->w[0][6] * lam2 / lam1 - 2 * (s->ndim - 2) * s->w[0][7];
  s->w[0][3] =
    (1 -
     ((lam0 + lam1 + lam2) / 7 -
      (lam0 * lam1 + lam0 * lam2 + lam1 * lam2) / 5) * 9 -
     lam0 * 3 * lam1 * lam2) / (lam3 * 18 * (lam3 - lam0) * (lam3 -
							     lam1) * (lam3 -
								      lam2));
  s->w[0][2] =
    (1 -
     ((lam0 + lam1 + lam3) / 7 -
      (lam0 * lam1 + lam0 * lam3 + lam1 * lam3) / 5) * 9 -
     lam0 * 3 * lam1 * lam3) / (lam2 * 18 * (lam2 - lam0) * (lam2 -
							     lam1) * (lam2 -
								      lam3)) -
    2 * (s->ndim - 1) * s->w[0][6];
  s->w[0][1] =
    (1 -
     ((lam0 + lam2 + lam3) / 7 -
      (lam0 * lam2 + lam0 * lam3 + lam2 * lam3) / 5) * 9 -
     lam0 * 3 * lam2 * lam3) / (lam1 * 18 * (lam1 - lam0) * (lam1 -
							     lam2) * (lam1 -
								      lam3)) -
    2 * (s->ndim - 1) * (s->w[0][6] + s->w[0][5] +
			 (s->ndim - 2) * s->w[0][7]);

/*     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules */
  s->w[1][s->wtleng - 1] = 1 / (108 * lam0 * lam0 * lam0 * lam0) / twondm;
  if (s->ndim > 2) {
    d = (6 * lam1) * (6 * lam1) * (6 * lam1);
    s->w[1][7] = (1 - 27 * twondm * s->w[1][8] * lam0 * lam0 * lam0) / d;
  }
  s->w[1][6] =
    (1 - lam1 * 5 / 3 -
     twondm * 15 * s->w[1][s->wtleng - 1] * (lam0 * lam0) * (lam0 -
							     lam1)) / (lam1 *
								       60 *
								       lam2 *
								       (lam2 -
									lam1));
  s->w[1][5] =
    (1 -
     9 * (8 * lam1 * lam2 * s->w[1][6] +
	  twondm * s->w[1][s->wtleng -
			   1] * lam0 * lam0)) / (lam1 * 36 * lam1) -
    2 * s->w[1][7] * (s->ndim - 2);
  s->w[1][3] =
    (1 -
     7 * ((lam1 + lam2) / 5 - lam1 * lam2 / 3 +
	  twondm * s->w[1][s->wtleng - 1] * lam0 * (lam0 - lam1) * (lam0 -
								    lam2))) /
    (14 * lam3 * (lam3 - lam1) * (lam3 - lam2));
  s->w[1][2] =
    (1 -
     7 * ((lam1 + lam3) / 5 - lam1 * lam3 / 3 +
	  twondm * s->w[1][s->wtleng - 1] * lam0 * (lam0 - lam1) * (lam0 -
								    lam3))) /
    (lam2 * 14 * (lam2 - lam1) * (lam2 - lam3)) - 2 * (s->ndim -
						       1) * s->w[1][6];
  s->w[1][1] =
    (1 -
     7 * ((lam2 + lam3) / 5 - lam2 * lam3 / 3 +
	  twondm * s->w[1][s->wtleng - 1] * lam0 * (lam0 - lam2) * (lam0 -
								    lam3))) /
    (lam1 * 14 * (lam1 - lam2) * (lam1 - lam3)) - 2 * (s->ndim -
						       1) * (s->w[1][6] +
							     s->w[1][5] +
							     (s->ndim -
							      2) *
							     s->w[1][7]);

  s->w[2][s->wtleng - 1] = 5 / (324 * lam0 * lam0 * lam0 * lam0) / twondm;
  if (s->ndim > 2)
    s->w[2][7] = (1 - 27 * twondm * s->w[2][8] * lam0 * lam0 * lam0) / d;

  s->w[2][6] =
    (1 - lam1 * 5 / 3 -
     twondm * 15 * s->w[2][s->wtleng - 1] * (lam0 * lam0) * (lam0 -
							     lam1)) / (lam1 *
								       60 *
								       lam2 *
								       (lam2 -
									lam1));
  s->w[2][5] =
    (1 -
     9 * (lam1 * 8 * lam2 * s->w[2][6] +
	  twondm * s->w[2][s->wtleng -
			   1] * lam0 * lam0)) / (lam1 * 36 * lam1) -
    s->w[2][7] * 2 * (s->ndim - 2);
  s->w[2][4] =
    (1 -
     7 * ((lam1 + lam2) / 5 - lam1 * lam2 / 3 +
	  twondm * s->w[2][s->wtleng - 1] * lam0 * (lam0 - lam1) * (lam0 -
								    lam2))) /
    (lamp * 14 * (lamp - lam1) * (lamp - lam2));
  s->w[2][2] =
    (1 -
     7 * ((lam1 + lamp) / 5 - lam1 * lamp / 3 +
	  twondm * s->w[2][s->wtleng - 1] * lam0 * (lam0 - lam1) * (lam0 -
								    lamp))) /
    (lam2 * 14 * (lam2 - lam1) * (lam2 - lamp)) - 2 * (s->ndim -
						       1) * s->w[2][6];
  s->w[2][1] =
    (1 -
     7 * ((lam2 + lamp) / 5 - lam2 * lamp / 3 +
	  twondm * s->w[2][s->wtleng - 1] * lam0 * (lam0 - lam2) * (lam0 -
								    lamp))) /
    (lam1 * 14 * (lam1 - lam2) * (lam1 - lamp)) - 2 * (s->ndim -
						       1) * (s->w[2][6] +
							     s->w[2][5] +
							     (s->ndim -
							      2) *
							     s->w[2][7]);

  s->w[3][s->wtleng - 1] = 2 / (81 * lam0 * lam0 * lam0 * lam0) / twondm;
  if (s->ndim > 2)
    s->w[3][7] = (2 - 27 * twondm * s->w[3][8] * lam0 * lam0 * lam0) / d;

  s->w[3][6] =
    (2 - lam1 * 15 / 9 -
     twondm * 15 * s->w[3][s->wtleng - 1] * lam0 * (lam0 -
						    lam1)) / (lam1 * 60 *
							      lam2 * (lam2 -
								      lam1));
  s->w[3][5] =
    (1 -
     9 * (lam1 * 8 * lam2 * s->w[3][6] +
	  twondm * s->w[3][s->wtleng -
			   1] * (lam0 * lam0))) / (lam1 * 36 * lam1) -
    s->w[3][7] * 2 * (s->ndim - 2);
  s->w[3][3] =
    (2 -
     7 * ((lam1 + lam2) / 5 - lam1 * lam2 / 3 +
	  twondm * s->w[3][s->wtleng - 1] * lam0 * (lam0 - lam1) * (lam0 -
								    lam2))) /
    (lam3 * 14 * (lam3 - lam1) * (lam3 - lam2));
  s->w[3][2] =
    (2 -
     7 * ((lam1 + lam3) / 5 - lam1 * lam3 / 3 +
	  twondm * s->w[3][s->wtleng - 1] * lam0 * (lam0 - lam1) * (lam0 -
								    lam3))) /
    (lam2 * 14 * (lam2 - lam1) * (lam2 - lam3)) - 2 * (s->ndim -
						       1) * s->w[3][6];
  s->w[3][1] =
    (2 -
     7 * ((lam2 + lam3) / 5 - lam2 * lam3 / 3 +
	  twondm * s->w[3][s->wtleng - 1] * lam0 * (lam0 - lam2) * (lam0 -
								    lam3))) /
    (lam1 * 14 * (lam1 - lam2) * (lam1 - lam3)) - 2 * (s->ndim -
						       1) * (s->w[3][6] +
							     s->w[3][5] +
							     (s->ndim -
							      2) *
							     s->w[3][7]);
  s->w[4][1] = 1 / (lam1 * 6);


/*     Set generator values */

  lam0 = sqrt (lam0);
  lam1 = sqrt (lam1);
  lam2 = sqrt (lam2);
  lam3 = sqrt (lam3);
  lamp = sqrt (lamp);

  // fprintf (stderr, "lam0 = %10.7f  %10.7f %10.7f %10.7f %10.7f \n", lam0,lam1,lam2,lam3,lamp);

  for (i = 0; i < s->ndim; i++)
    s->g[i][s->wtleng - 1] = lam0;

  if (s->ndim > 2) {
    s->g[0][7] = lam1;
    s->g[1][7] = lam1;
    s->g[2][7] = lam1;
  }
  s->g[0][6] = lam1;
  s->g[1][6] = lam2;
  s->g[0][5] = lam1;
  s->g[1][5] = lam1;
  s->g[0][4] = lamp;
  s->g[0][3] = lam3;
  s->g[0][2] = lam2;
  s->g[0][1] = lam1;


/*     Compute final weight values. */

/*     The null rule weights are computed from differences between */

/*     the degree 9 rule weights and lower degree rule weights. */

  s->w[0][0] = twondm;
  for (j = 1; j < 5; j++) {
    for (i = 1; i < s->wtleng; i++) {
      s->w[j][i] = s->w[j][i] - s->w[0][i];
      s->w[j][0] = s->w[j][0] - s->rulpts[i] * s->w[j][i];
    }
  }

  for (i = 1; i < s->wtleng; i++) {
    s->w[0][i] = twondm * s->w[0][i];
    s->w[0][0] -= s->rulpts[i] * s->w[0][i];
  }

  s->errcof[0] = 5.;
  s->errcof[1] = 5.;
  s->errcof[2] = 1.f;
  s->errcof[3] = 5.;
  s->errcof[4] = .5f;
  s->errcof[5] = .25f;

  return 0;
}				/* End of  d09hre_ */



int
d07hre_ (dcuhre_state * s)
{


  /* Local variables */
  int i, j;
  double lam0, lam1, lam2, lamp, ratio, twondm;

  for (j = 0; j < s->wtleng; j++) {
    for (i = 0; i < s->ndim; i++) {
      s->w[i][j] = 0.0;
    }
    for (i = 0; i < 5; i++) {
      s->w[i][j] = 0.0;
    }
    s->rulpts[j] = 2 * s->ndim;
  }

  twondm = (double) pow_i (2, s->ndim);
  s->rulpts[s->wtleng - 1] = twondm;
  s->rulpts[s->wtleng - 2] = (double) (2 * s->ndim * (s->ndim - 1));
  s->rulpts[0] = 1.;

/*     Compute squared generator parameters */

  lam0 = .4707f;
  lamp = .5625f;
  lam1 = 4 / (15 - 5 / lam0);
  ratio = (1 - lam1 / lam0) / 27;
  lam2 =
    (5 - lam1 * 7 - ratio * 35) / (7 - lam1 * 35 / 3 - ratio * 35 / lam0);

/*     Compute degree 7 rule weights */

  s->w[0][5] = 1 / (3 * lam0) / (3 * lam0) / (3 * lam0) / twondm;
  s->w[0][4] = (1 - lam0 * 5 / 3) / ((lam1 - lam0) * 60 * lam1 * lam1);
  s->w[0][2] =
    (1 - lam2 * 5 / 3 -
     twondm * 5 * s->w[0][5] * lam0 * (lam0 - lam2)) / (lam1 * 10 * (lam1 -
								     lam2)) -
    2 * (s->ndim - 1) * s->w[0][4];
  s->w[0][1] =
    (1 - lam1 * 5 / 3 -
     twondm * 5 * s->w[0][5] * lam0 * (lam0 - lam1)) / (lam2 * 10 * (lam2 -
								     lam1));

/*     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules */

  s->w[1][5] = 1 / (lam0 * lam0 * lam0 * 36) / twondm;
  s->w[1][4] =
    (1 - twondm * 9 * s->w[1][5] * lam0 * lam0) / (lam1 * lam1 * 36);
  s->w[1][2] =
    (1 - lam2 * 5 / 3 -
     twondm * 5 * s->w[1][5] * lam0 * (lam0 - lam2)) / (lam1 * 10 * (lam1 -
								     lam2)) -
    2 * (s->ndim - 1) * s->w[1][4];
  s->w[1][1] =
    (1 - lam1 * 5 / 3 -
     twondm * 5 * s->w[1][5] * lam0 * (lam0 - lam1)) / (lam2 * 10 * (lam2 -
								     lam1));
  s->w[2][5] = 5 / (lam0 * lam0 * lam0 * 108) / twondm;
  s->w[2][4] =
    (1 - twondm * 9 * s->w[2][5] * lam0 * lam0) / (lam1 * lam1 * 36);
  s->w[2][2] =
    (1 - lamp * 5 / 3 -
     twondm * 5 * s->w[2][5] * lam0 * (lam0 - lamp)) / (lam1 * 10 * (lam1 -
								     lamp)) -
    2 * (s->ndim - 1) * s->w[2][4];
  s->w[2][3] =
    (1 - lam1 * 5 / 3 -
     twondm * 5 * s->w[2][5] * lam0 * (lam0 - lam1)) / (lamp * 10 * (lamp -
								     lam1));
  s->w[3][5] = 1 / (lam0 * lam0 * lam0 * 54) / twondm;
  s->w[3][4] =
    (1 - twondm * 18 * s->w[3][5] * lam0 * lam0) / (lam1 * lam1 * 72);
  s->w[3][2] =
    (1 - lam2 * 10 / 3 -
     twondm * 10 * s->w[3][5] * lam0 * (lam0 - lam2)) / (lam1 * 20 * (lam1 -
								      lam2)) -
    2 * (s->ndim - 1) * s->w[3][4];
  s->w[3][1] =
    (1 - lam1 * 10 / 3 -
     twondm * 10 * s->w[3][5] * lam0 * (lam0 - lam1)) / (lam2 * 20 * (lam2 -
								      lam1));



/*     Set generator values */

  lam0 = sqrt (lam0);
  lam1 = sqrt (lam1);
  lam2 = sqrt (lam2);
  lamp = sqrt (lamp);
  for (i = 0; i < s->ndim; i++) {
    s->g[i][s->wtleng - 1] = lam0;
  }
  s->g[0][s->wtleng - 2] = lam1;
  s->g[1][s->wtleng - 2] = lam1;
  s->g[0][s->wtleng - 5] = lam2;
  s->g[0][s->wtleng - 4] = lam1;
  s->g[0][s->wtleng - 3] = lamp;

/*     Compute final weight values. */

/*     The null rule weights are computed from differences between */

/*     the degree 7 rule weights and lower degree rule weights. */

  s->w[0][0] = twondm;
  for (j = 1; j < 5; j++) {
    for (i = 1; i < s->wtleng; i++) {
      s->w[j][i] -= s->w[0][i];
      s->w[j][0] -= s->rulpts[i] * s->w[j][i];
    }
  }

  for (i = 1; i < s->wtleng; i++) {
    s->w[0][i] = twondm * s->w[0][i];
    s->w[0][0] -= s->rulpts[i] * s->w[0][i];
  }

/*     Set error coefficients */
  s->errcof[0] = 5.;
  s->errcof[1] = 5.;
  s->errcof[2] = 1.;
  s->errcof[3] = 5.;
  s->errcof[4] = .5f;
  s->errcof[5] = .25f;


  return 0;
}				/* End of d07hre_ */









int
dchhre_ (dcuhre_state * s)
{
  int j;

  /*   Check on legal KEY. */

  if (s->key < 0 || s->key > 4) {
    s->ifail = 2;
    return 1;
  }

  /*   Check on legal NDIM. */
  if (s->ndim < 2 || s->ndim > maxdim) {
    s->ifail = 3;
    return 1;
  }

  /*   For KEY = 1, NDIM must be equal to 2. */
  if (s->key == 1 && s->ndim != 2) {
    s->ifail = 4;
    return 1;
  }

  /*   For KEY = 2, NDIM must be equal to 3. */
  if (s->key == 2 && s->ndim != 3) {
    s->ifail = 5;
    return 1;
  }

  /*   For KEY = 0, we point at the selected integration rule. */
  if (s->key == 0) {
    if (s->ndim == 2) {
      s->keyf = 1;
    } else if (s->ndim == 3) {
      s->keyf = 2;
    } else {
      s->keyf = 3;
    }
  } else {
    s->keyf = s->key;
  }

  /*   Compute NUM and WTLENG as a function of KEYF and NDIM. */
  if (s->keyf == 1) {
    s->num = 65;
    s->wtleng = 14;
  } else if (s->keyf == 2) {
    s->num = 127;
    s->wtleng = 13;
  } else if (s->keyf == 3) {
    s->num =
      1 + 4 * 2 * s->ndim + 2 * s->ndim * (s->ndim - 1) +
      4 * s->ndim * (s->ndim - 1) + 4 * s->ndim * (s->ndim - 1) * (s->ndim -
								   2) / 3 +
      pow_i (2, s->ndim);
    s->wtleng = 9;
    if (s->ndim == 2) {
      s->wtleng = 8;
    }
  } else if (s->keyf == 4) {
    s->num =
      1 + s->ndim * 6 + (s->ndim * 2) * (s->ndim - 1) + pow_i (2, s->ndim);
    s->wtleng = 6;
  }



  /*   Compute MAXSUB. */
  //s->maxsub = (s->maxcls - s->num) / (s->num *2) + 1;
  s->maxsub = (s->maxcls - s->num) / (s->num) + 1;
  /*   Compute MINSUB. */
  s->minsub = (s->mincls - s->num) / (s->num * 2) + 1;
  if ((s->mincls - s->num) % (s->num << 1) != 0) {
    s->minsub++;
  }
  s->minsub = max (2, s->minsub);


  /*   Check on legal upper and lower limits of integration. */
  for (j = 0; j < s->ndim; ++j) {
    if (s->xl[j] - s->xu[j] >= 0.) {
      fprintf (stderr, "%d  %f    %f\n", j, s->xl[j], s->xu[j]);
      s->ifail = 7;
      return 0;
    }
  }

  /*   Check on MAXPTS < 3*NUM. */
  if (s->maxcls < s->num * 3) {
      fprintf(stderr,"Increase maxcls =%d to ", s->maxcls);
    s->maxcls += s->num * 3; 
    //    fprintf(stderr,"ifail =8   with maxcls =%d  num =%d \n",s->maxcls, s->num );
      fprintf(stderr,"New maxcls =%d\n", s->maxcls);
      //    s->ifail = 8;
      //    return 0;
  }

  /*   Check on MAXPTS >= MINPTS. */
  if (s->maxcls < s->mincls) {
    s->ifail = 9;
    return 1;
  }

  /*   Check on legal accuracy requests. */
  if (s->epsabs < 0. && s->epsrel < 0.) {
    s->ifail = 10;
    return 0;
  }

  /*    Check on legal RESTAR. */
  if (s->restar != 0 && s->restar != 1) {
    s->ifail = 12;
    return 1;
  }

  return 0;

/* ***END DCHHRE */
}



int
initialize_state (dcuhre_state * s, double *a, double *b, int dim)
{
  /* This function initializes part of the global variables which can be initialized 
     without any other information.
     More dependent variables are initialized in dcuhre_ function after dchhre_
   */

  s->ndim = dim;
  s->xl = a;
  s->xu = b;


  /* default */
  s->mincls = 0;
  s->maxcls = 10000;
  s->key = 0;
  s->epsabs = 0.0;
  s->epsrel = 0.01;

  /*if (dim==4){
     s->maxcls=2000;
     s->epsrel=0.01;
     }else if (dim>4){
     s->maxcls=10000;
     s->epsrel=0.01;  
     }else{
     s->maxcls=1000;
     s->epsrel=0.01;  
     } */

  s->result = 0.0;
  s->error = 0.0;
  s->total_neval = 0;
  s->ifail = 0;

  s->restar = 0;
  s->verbose = 0;
  s->sbrgns = 0;

  s->numfun = 1;		/*change this for vector functions */

  s->vol_rate =1.0;

  s->scale=0;

  return 0;
}


/* integer pow function  */
int
pow_i (int base, int expo)
{

  int i, res = 1;

  for (i = 0; i < expo; i++) {
    res *= base;
  }
  return res;
}



int
sbrg_free (sub_region ** sbrg_heap, int maxsub)
{

  int i = 0;
  sub_region *cw_sbrg;

  while (sbrg_heap[i] != NULL && i < maxsub) {
    cw_sbrg = sbrg_heap[i];
    free (cw_sbrg->center);
    free (cw_sbrg->hwidth);
    free (sbrg_heap[i++]);
  }

  return 0;
}


int
print_rule (dcuhre_state * s)
{

  int i, j;

  fprintf (stderr, "Number of function evaluation for each subregion is %d\n",
	   s->num);
  fprintf (stderr, "Printing g whose size is dim=%d x wtleng =%d\n", s->ndim,
	   s->wtleng);
  for (i = 0; i < s->ndim; i++) {
    for (j = 0; j < s->wtleng; j++)
      fprintf (stderr, " %f", s->g[i][j]);
    fprintf (stderr, "\n");
  }

  fprintf (stderr, "Printing w whose size is 5 x wtleng =%d\n", s->wtleng);
  for (i = 0; i < 5; i++) {
    for (j = 0; j < s->wtleng; j++)
      fprintf (stderr, " %f", s->w[i][j]);
    fprintf (stderr, "\n");
  }

  fprintf (stderr, "Printing rulpts whose size is  wtleng =%d\n", s->wtleng);
  for (j = 0; j < s->wtleng; j++)
    fprintf (stderr, " %d", s->rulpts[j]);
  fprintf (stderr, "\n");


  fprintf (stderr, "Printing scales whose size is 3 x wtleng =%d\n",
	   s->wtleng);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < s->wtleng; j++)
      fprintf (stderr, " %f", s->scales[i][j]);
    fprintf (stderr, "\n");
  }

  fprintf (stderr, "Printing norms whose size is s x wtleng =%d\n",
	   s->wtleng);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < s->wtleng; j++)
      fprintf (stderr, " %f", s->norms[i][j]);
    fprintf (stderr, "\n");
  }


  fprintf (stderr, "Printing rulpts whose size is 6\n");
  for (j = 0; j < 6; j++)
    fprintf (stderr, " %f", s->errcof[j]);
  fprintf (stderr, "\n");

  return 0;

}


int
print_sbrg (sub_region * cw_sbrg, int dim)
{
  int i;

  fprintf (stderr, "Printing subregion with id = %d\n", cw_sbrg->region_id);
  fprintf (stderr, "                     level = %d\n",
	   cw_sbrg->region_level);
  fprintf (stderr, "                left child = %d right child =%d\n",
	   cw_sbrg->lchild_id, cw_sbrg->rchild_id);
  fprintf (stderr, "               local error = %f\n", cw_sbrg->local_error);
  fprintf (stderr, "              local result = %f\n",cw_sbrg->local_result);
  fprintf (stderr, "                     scale = %f\n",cw_sbrg->cur_scale);
  fprintf (stderr, "on the sub region center   hwidth\n");
  for (i = 0; i < dim; i++) {
    fprintf (stderr, "                 %f, %f\n", cw_sbrg->center[i],
	     cw_sbrg->hwidth[i]);
  }

  return 0;
}
