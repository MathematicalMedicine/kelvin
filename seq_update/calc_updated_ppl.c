/*
 * calc_updated_ppl - Calculate PPL, LD-PPL and PPLD from Kelvin-format
 *                    sequentially updated 2-point average Het LR files.
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright 2007, The Research Institute at Nationwide Children's Hospital
 * All rights reserved. Permission is granted to use this software for
 * non-profit educational purposes only.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>

/* TO DO:
 * - sex-specific theta cutoffs
 * - seperate linkage prior / ld prior
 */

#define DPRIME_COL 1
#define THETA_COL  2
#define LR_COL     3

typedef struct {
  char name1[32],
    name2[32];
  int num,
    numcols,
    numdprimes,
    numthetas,
    no_ld,
    holey_grid,
    *datacols;
} st_marker;

typedef struct {
  float *dprimes,
    *thetas;
  double lr;
} st_data;

typedef struct {
  int arrsize,
    numelems,
    lastidx,
    dimsize;
  float *arr;
} st_dim;

typedef struct {
  short numdims;
  int totalelems;
  st_dim *dims;
} st_multidim;


/* Default values that can be overridden on the command line */
int sexspecific = 0;   /* -s sex-specific thetas */
double prior = 0.02;   /* -p prior probability */
double weight = 0.95;  /* -w weighting factor */
double cutoff = 0.05;  /* -c theta cutoff */

/* Globally available for error messages */
int lineno = 0;

double calc_ppl_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr);
double calc_ppl_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr);
double calc_ldppl_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr);
double calc_ldppl_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr);
double calc_ppld_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr);
double calc_ppld_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr);
int parse_command_line (int argc, char **argv);
int get_marker_line (st_marker *marker, FILE *fp);
int get_header_line (st_marker *marker, st_data *data, FILE *fp);
int get_data_line (st_marker *marker, st_data *data, FILE *fp);
int multi_insert (st_multidim *md, float *vals, int num);
int multi_find (st_multidim *md, float *vals, int num);
int insert (st_dim *dim, float val);
int find (st_dim *dim, float val, int *found);
void multi_dump (st_multidim *md);

int main (int argc, char **argv)
{
  int argidx, ret, datalines, va, vb;
  FILE *fp;
  st_marker marker;
  st_data data;
  st_multidim dprimes, thetas;
  double **lr;

  argidx = parse_command_line (argc, argv);
  if (argidx >= argc) {
    fprintf (stderr, "missing file name\n");
    exit (-1);
  }
  if ((fp = fopen (argv[argidx], "r")) == NULL) {
    fprintf (stderr, "open '%s' failed, %s\n", argv[argidx], strerror (errno));
    exit (-1);
  }
  
  memset (&marker, 0, sizeof (st_marker));
  memset (&dprimes, 0, sizeof (st_multidim));
  memset (&thetas, 0, sizeof (st_multidim));

  if ((ret = get_marker_line (&marker, fp)) == 0) {
    fprintf (stderr, "file '%s' is empty\n", argv[argidx]);
    exit (-1);
  } else if (ret == -1) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", lineno, argv[argidx]);
    exit (-1);
  }
  if ((ret = get_header_line (&marker, &data, fp)) == 0) {
    fprintf (stderr, "file '%s' ends unexpectedly at line %d\n", argv[argidx], lineno);
    exit (-1);
  } else if (ret == -1) {
    fprintf (stderr, "can't parse header, line %d in file '%s'\n", lineno, argv[argidx]);
    exit (-1);
  }

  while ((ret = get_data_line (&marker, &data, fp)) == 1) {
    if (multi_insert (&dprimes, data.dprimes, marker.numdprimes) == -1) {
      fprintf (stderr, "insert into dprimes failed, %s\n", strerror (errno));
      exit (-1);
    }
    if (multi_insert (&thetas, data.thetas, marker.numthetas) == -1) {
      fprintf (stderr, "insert into thetas failed, %s\n", strerror (errno));
      exit (-1);
    }
  }
  if (ret == -1) {
    fprintf (stderr, "can't parse data, line %d in file '%s'\n", lineno, argv[argidx]);
    exit (-1);
  }

  if (fseek (fp, 0, SEEK_SET) == -1) {
    fprintf (stderr, "fseek failed, %s\n", strerror (errno));
    exit (-1);
  }
  lineno = 0;

  if ((lr = malloc (sizeof (double *) * dprimes.totalelems)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes.totalelems; va++) {
    if ((lr[va] = malloc (sizeof (double) * thetas.totalelems)) == NULL) {
      fprintf (stderr, "malloc failed, %s\n", strerror (errno));
      exit (-1);
    }
  }

  if (! marker.no_ld)
    printf ("%34s %6s %6s %6s\n", " ", "PPL","LD-PPL", "PPLD");
  else
    printf ("%34s %6s\n", " ", "PPL");

  while ((ret = get_marker_line (&marker, fp)) != 0) {
    if ((ret = get_header_line (&marker, &data, fp)) == 0) {
      fprintf (stderr, "file '%s' ends unexpectedly at line %d\n", argv[argidx], lineno);
      exit (-1);
    } else if (ret == -1) {
      fprintf (stderr, "can't parse header, line %d in file '%s'\n", lineno, argv[argidx]);
      exit (-1);
    }
    datalines = 0;
    for (va = 0; va < dprimes.totalelems; va++) 
      for (vb = 0; vb < thetas.totalelems; vb++)
	lr[va][vb] = DBL_MAX;
    while ((ret = get_data_line (&marker, &data, fp)) == 1) {

      if ((va = multi_find (&dprimes, data.dprimes, marker.numdprimes)) == -1) {
	fprintf (stderr, "unexpected dprimes at line %d\n", lineno);
	exit (-1);
      }

      if ((vb = multi_find (&thetas, data.thetas, marker.numthetas)) == -1) {
	fprintf (stderr, "unexpected thetas at line %d\n", lineno);
	exit (-1);
      }
      if (lr[va][vb] != DBL_MAX) {
	fprintf (stderr, "duplicate dprime/theta at line %d\n", lineno);
	fprintf (stderr, "orig value is %.4f (%d/%d)\n", lr[va][vb], va, vb);
	exit (-1);
      }
      lr[va][vb] = data.lr;
      datalines++;
    }
    if (ret == -1) {
      fprintf (stderr, "can't parse data, line %d in file '%s'\n", lineno, argv[argidx]);
      exit (-1);
    }
    if ((! marker.holey_grid) && (datalines != (va = dprimes.totalelems * thetas.totalelems))) {
      fprintf (stderr, "expected %d data lines, found %d, ending at line %d\n", va, 
	       datalines, lineno);
      exit (-1);
    }
    printf ("%2d %15s %15s", marker.num, marker.name1, marker.name2);
    printf (" %6.4f", (! sexspecific) ? calc_ppl_sexavg (&dprimes, &thetas, lr) : 
	    calc_ppl_sexspc (&dprimes, &thetas, lr));
    if (! marker.no_ld) {
      printf (" %6.4f", (! sexspecific) ? calc_ldppl_sexavg (&dprimes, &thetas, lr) : 
	      calc_ldppl_sexspc (&dprimes, &thetas, lr));
      printf (" %6.4f", (! sexspecific) ? calc_ppld_sexavg (&dprimes, &thetas, lr) : 
	      calc_ppld_sexspc (&dprimes, &thetas, lr));
    }
    printf ("\n");
  }
  if (ret == -1) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", lineno, argv[argidx]);
    exit (-1);
  }
  
  for (va = 0; va < dprimes.totalelems; va++) {
    free (lr[va]);
  }
  free (lr);

  for (va = 0; va < dprimes.numdims; va++)
    free (dprimes.dims[va].arr);
  free (dprimes.dims);
  for (va = 0; va < thetas.numdims; va++)
    free (thetas.dims[va].arr);
  free (thetas.dims);

  free (data.dprimes);
  free (data.thetas);
  free (marker.datacols);

  exit (0);
}

/* For sex-averaged thetas, we calculate the area under a two-dimensional
 * curve defined by theta (X axis) vs. LR (Y axis), where D' (all D' if
 * more than one) is 0. The area is calculated by computing the areas of
 * polygons bounded by adjascent thetas and associated LRs.
 */
double calc_ppl_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr)
{
  int va, vb;
  float *zeros;
  double lr1, lr2, mtheta1, mtheta2, cutlr, pre_int=0, post_int=0, integral, ppl;
  st_dim *mthetas;
  
  if ((zeros = malloc (sizeof (float) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((va = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);

  /* if any LR where D' == 0 is FLT_MAX (that is, undefined), then there
   * shouldn't be any defined LR where D' == 0, regardless of theta, and
   * we can't calculate a PPL.
   */
  if (lr[va][0] == FLT_MAX)
    return (-1);

  mthetas = &thetas->dims[0];
  for (vb = 1; vb < mthetas->numelems; vb++) {
    mtheta1 = mthetas->arr[vb-1];
    mtheta2 = mthetas->arr[vb];
    lr1 = lr[va][vb-1];
    lr2 = lr[va][vb];

    if (mtheta1 >= cutoff) {
      /* entire region outside the cutoff */
      post_int += (lr1 + lr2) * (mtheta2 - mtheta1);
      
    } else if (mtheta2 <= cutoff) {
      /* entire region inside the cutoff */
      pre_int += (lr1 + lr2) * (mtheta2 - mtheta1);
      
    } else {
      /* region straddles the cutoff */
      cutlr = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
      pre_int += (lr1 + cutlr) * (cutoff - mtheta1);
      post_int += (cutlr + lr2) * (mtheta2 - cutoff);
    }
  }

  /* The 0.5 here is factored out of the calculation of area. */
  pre_int *= (weight / cutoff) * 0.5;
  post_int *= ((1 - weight) / (0.5 - cutoff)) * 0.5;
  integral = pre_int + post_int;
  ppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ppl);
}


/* For sex-specific thetas, we calculate the volume under a three-dimensional
 * curve defined by male theta (X axis) vs. female theta (Y axis) vs. LR (Z
 * axis), where D' (all D', if more than one) is 0. The volume is calculated
 * by computing the volumes of polyhedrons bounded by adjascent males thetas,
 * adjascent female thetas, and the four LRs associated with the intersections
 * of the male and female thetas. Looking at the two-dimensional X-Y plane,
 * LR1 would be in the lower left, LR2 in the lower right, LR3 in the upper
 * left, and LR4 in the upper right.
 */
double calc_ppl_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr)
{
  int va, vb, vc;
  float *zeros;
  double lr1, lr2, lr3, lr4, mtheta1, mtheta2, ftheta1, ftheta2, cutlr2, cutlr3, cutlr4;
  double pre_int=0, post_int=0, cutvol, integral, ppl;
  st_dim *mthetas, *fthetas;
  
  if ((zeros = malloc (sizeof (float) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((va = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);

  /* See comment in calc_ppl_sexavg() */
  if (lr[va][0] == FLT_MAX)
    return (-1);
  
  /* Does it really matter which one is male or female? */
  mthetas = &thetas->dims[0];
  fthetas = &thetas->dims[1];
  for (vb = 1 ; vb < mthetas->numelems; vb++) {
    mtheta1 = mthetas->arr[vb - 1];
    mtheta2 = mthetas->arr[vb];
    for (vc = 1 ; vc < fthetas->numelems; vc++) {
      ftheta1 = fthetas->arr[vc - 1];
      ftheta2 = fthetas->arr[vc];
      lr1 = lr[va][(vb - 1) * mthetas->numelems + vc - 1];
      lr2 = lr[va][vb * mthetas->numelems + vc - 1];
      lr3 = lr[va][(vb - 1) * mthetas->numelems + vc];
      lr4 = lr[va][vb * mthetas->numelems + vc];
      
      if ((mtheta1 >= cutoff) || (ftheta1 >= cutoff)) {
	/* entire region is outside cutoff area */
	post_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	
      } else if ((mtheta2 <= cutoff) && (ftheta2 <= cutoff)) {
	/* entire region inside cutoff area */
	pre_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	
      } else if ((ftheta2 <= cutoff) && ((mtheta1 <= cutoff) && (cutoff <= mtheta2))) {
	/* region straddles cutoff on male axis */
	cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	pre_int += (cutoff - mtheta1) * (ftheta2 - ftheta1) * (lr1 + cutlr2 + lr3 + cutlr4);
	post_int += (mtheta2 - cutoff) * (ftheta2 - ftheta1) * (cutlr2 + lr2 + cutlr4 + lr4);
	
      } else if ((mtheta2 <= cutoff) && ((ftheta1 <= cutoff) && (cutoff <= ftheta2))) {
	/* region straddles cutoff on female axis */
	cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	cutlr4 = lr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr4 - lr2);
	pre_int += (cutoff - ftheta1) * (mtheta2 - mtheta1) * (lr1 + lr2 + cutlr3 + cutlr4);
	post_int += (ftheta2 - cutoff) * (mtheta2 - mtheta1) * (cutlr3 + cutlr4 + lr3 + lr4);
	
      } else {
	/* region straddles cutoff on both axes */
	post_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	cutlr4 = cutlr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (cutlr4 - cutlr2);
	cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	cutvol = (cutoff - mtheta1) * (cutoff - ftheta1) * (lr1 + cutlr2 + cutlr3 + cutlr4);
	post_int -= cutvol;
	pre_int += cutvol;
      }
    }
  }

  /* The 0.25 here is factored out of the caclculation of volume */
  pre_int *= (weight / (cutoff * cutoff)) * 0.25;
  post_int *= ((1 - weight) / (0.25 - (cutoff * cutoff))) * 0.25;
  integral = pre_int + post_int;
  ppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ppl);
}


/* This is just like calc_ppl_sexavg, except that instead of only 
 * calculating with the LR where D' is 0, we average the LR across
 * all D', and calculate the area under the curve defined by 
 * theta vs. avg LR.
 */
double calc_ldppl_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr)
{
  int va, vb, numdprimes=0;
  double lr1, lr2=0, mtheta1, mtheta2, cutlr, pre_int=0, post_int=0, integral, ldppl;
  st_dim *mthetas;

  mthetas = &thetas->dims[0];
  for (va = 0; va < dprimes->totalelems; va++) {
    if (lr[va][0] == FLT_MAX)
      continue;
    lr2 += lr[va][0];
    numdprimes++;
  }
  lr2 /= (double) numdprimes;

  for (vb = 1; vb < mthetas->numelems; vb++) {
    mtheta1 = mthetas->arr[vb-1];
    mtheta2 = mthetas->arr[vb];
    lr1 = lr2;
    lr2 = 0;
    for (va = 0; va < dprimes->totalelems; va++) {
      if (lr[va][vb] == FLT_MAX)
	continue;
      lr2 += lr[va][vb];
    }
    lr2 /= (double) numdprimes;
    
    if (mtheta1 >= cutoff) {
      post_int += (lr1 + lr2) * (mtheta2 - mtheta1);
      
    } else if (mtheta2 <= cutoff) {
      pre_int += (lr1 + lr2) * (mtheta2 - mtheta1);
      
    } else {
      cutlr = lr1 +
	((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
      pre_int += (lr1 + cutlr) * (cutoff - mtheta1);
      post_int += (cutlr + lr2) * (mtheta2 - cutoff);
    }
  }

  /* The 0.5 here is factored out of the calculation of area. */
  pre_int *= (weight / cutoff) * 0.5;
  post_int *= ((1 - weight) / (0.5 - cutoff)) * 0.5;
  integral = pre_int + post_int;
  ldppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ldppl);
}


/* This is just like calc_ppl_sexspc, except instead of calculating
 * the volume using the LRs where D' is 0, we average each LR across
 * all D', and calculate the volume under the cur defined by male 
 * theta vs. female theta vs. avg LR.
 */
double calc_ldppl_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr)
{
  int va, vb, vc, numdprimes=0;
  double lr1, lr2, lr3, lr4, mtheta1, mtheta2, ftheta1, ftheta2, cutlr2, cutlr3, cutlr4;
  double pre_int=0, post_int=0, cutvol, integral, ldppl, *avglrs;
  st_dim *mthetas, *fthetas;

  mthetas = &thetas->dims[0];
  fthetas = &thetas->dims[1];
  if ((avglrs = malloc (sizeof (double) * fthetas->numelems)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }

  /* Calculate the first element of avglrs outside the loop to establish
   * the value of numdprimes.
   */
  avglrs[0] = 0;
  for (va = 0; va < dprimes->totalelems; va++) {
    if (lr[va][0] == FLT_MAX)
      continue;
    avglrs[0] += lr[va][0];
    numdprimes++;
  }
  avglrs[0] /= (double) numdprimes;

  /* Now calculate the rest of avglrs
   */
  for (vc = 1; vc < fthetas->numelems; vc++) {
    avglrs[vc] = 0;
    for (va = 0; va < dprimes->totalelems; va++) 
      if (lr[va][vc] == FLT_MAX)
	continue;
      avglrs[vc] += lr[va][vc];
    avglrs[vc] /= (double) numdprimes;
  }

  for (vb = 1; vb < mthetas->numelems; vb++) {
    mtheta1 = mthetas->arr[vb - 1];
    mtheta2 = mthetas->arr[vb];
    lr4 = 0;
    for (va = 0; va < dprimes->totalelems; va++) {
      if (lr[va][vb * mthetas->numelems] == FLT_MAX)
	continue;
      lr4 += lr[va][vb * mthetas->numelems];
    }
    lr4 /= (double) numdprimes;
    
    for (vc = 1; vc < fthetas->numelems; vc++) {
      ftheta1 = fthetas->arr[vc - 1];
      ftheta2 = fthetas->arr[vc];
      lr1 = avglrs[vc - 1];
      lr2 = lr4;
      lr3 = avglrs[vc];
      lr4 = 0;
      for (va = 0; va < dprimes->totalelems; va++) {
	if (lr[va][vb * mthetas->numelems + vc] == FLT_MAX)
	  continue;
	lr4 += lr[va][vb * mthetas->numelems + vc];
      }
      lr4 /= (double) numdprimes;

      if ((mtheta1 >= cutoff) || (ftheta1 >= cutoff)) {
	/* entire region is outside cutoff area */
	post_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	
      } else if ((mtheta2 <= cutoff) && (ftheta2 <= cutoff)) {
	/* entire region inside cutoff area */
	pre_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	
      } else if ((ftheta2 <= cutoff) && ((mtheta1 <= cutoff) && (cutoff <= mtheta2))) {
	/* region straddles cutoff on male axis */
	cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	pre_int += (cutoff - mtheta1) * (ftheta2 - ftheta1) * (lr1 + cutlr2 + lr3 + cutlr4);
	post_int += (mtheta2 - cutoff) * (ftheta2 - ftheta1) * (cutlr2 + lr2 + cutlr4 + lr4);
	
      } else if ((mtheta2 <= cutoff) && ((ftheta1 <= cutoff) && (cutoff <= ftheta2))) {
	/* region straddles cutoff on female axis */
	cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	cutlr4 = lr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr4 - lr2);
	pre_int += (cutoff - ftheta1) * (mtheta2 - mtheta1) * (lr1 + lr2 + cutlr3 + cutlr4);
	post_int += (ftheta2 - cutoff) * (mtheta2 - mtheta1) * (cutlr3 + cutlr4 + lr3 + lr4);
	
      } else {
	/* region straddles cutoff on both axes */
	post_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	cutlr4 = cutlr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (cutlr4 - cutlr2);
	cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	cutvol = (cutoff - mtheta1) * (cutoff - ftheta1) * (lr1 + cutlr2 + cutlr3 + cutlr4);
	post_int -= cutvol;
	pre_int += cutvol;
      }

      avglrs[vc - 1] = lr2;
    }
    avglrs[vc - 1] = lr4;
  }
  free (avglrs);

  /* The 0.25 here is factored out of the caclculation of volume */
  pre_int *= (weight / (cutoff * cutoff)) * 0.25;
  post_int *= ((1 - weight) / (0.25 - (cutoff * cutoff))) * 0.25;
  integral = pre_int + post_int;
  ldppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ldppl);
}


double calc_ppld_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr)
{
  int va, vb, ndprimes = 0, zero_idx;
  float *zeros;
  double mtheta1, mtheta2, lr1, lr2, cutlr;
  double ld_pre=0, ld_post=0, le_pre=0, le_post=0, denom, ppld;
  st_dim *mthetas;
 
  if ((zeros = malloc (sizeof (float) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((zero_idx = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);
  mthetas = &thetas->dims[0];

  for (va = 0; va < dprimes->totalelems; va++) {

    /* See comment in calc_ppl_sexavg() */
    if (lr[va][0] == FLT_MAX)
      continue;

    if (va == zero_idx) {
      /* D' == 0 */
      for (vb = 1; vb < mthetas->numelems; vb++) {
	mtheta1 = mthetas->arr[vb - 1];
	mtheta2 = mthetas->arr[vb];
	lr1 = lr[va][vb-1];
	lr2 = lr[va][vb];
	
	if (mtheta1 >= cutoff) {
	  le_post += (lr1 + lr2) * (mtheta2 - mtheta1);
	} else if (mtheta2 <= cutoff) {
	  le_pre += (lr1 + lr2) * (mtheta2 - mtheta1);
	} else {
	  cutlr = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	  le_pre += (lr1 + cutlr) * (cutoff - mtheta1);
	  le_post += (cutlr + lr2) * (mtheta2 - cutoff);
	}
      }
    } else {
      /* D' != 0 */
      ndprimes++;
      for (vb = 1; vb < mthetas->numelems; vb++) {	
	mtheta1 = mthetas->arr[vb - 1];
	mtheta2 = mthetas->arr[vb];
	lr1 = lr[va][vb-1];
	lr2 = lr[va][vb];
	
	if (mtheta1 >= cutoff) {
	  /* ld_post += (lr1 + lr2) * (mtheta2 - mtheta1) */;
	} else if (mtheta2 <= cutoff) {
	  ld_pre += (lr1 + lr2) * (mtheta2 - mtheta1);
	} else {
	  cutlr = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	  ld_pre += (lr1 + cutlr) * (cutoff - mtheta1);
	  /*ld_post += (cutlr + lr2) * (mtheta2 - cutoff);*/
	}
      }
    }
  }

  /* The 0.5 here is factored out of the calculation of area */
  /* We may be mis-using prior here by making no distiction between linkage prior and ld prior */
  ld_pre *= prior * (weight / cutoff) * 0.5;
  ld_pre /= ndprimes;
  /*ld_post *= ((1 - weight) / (0.5 - cutoff)) * 0.5;*/
  /*ld_post /= ndprimes;*/
  le_pre *= (1 - prior) * (weight / cutoff) * 0.5;
  le_post *= ((1 - weight) / (0.5 - cutoff)) * 0.5;
  if ((denom = ld_pre + ld_post + le_pre + le_post) == 0.0)
    return (0.0);
  ppld = (ld_pre + ld_post) / denom;
  return (ppld);
}


double calc_ppld_sexspc (st_multidim *dprimes, st_multidim *thetas, double**lr)
{
  int va, vb, vc, ndprimes = 0, zero_idx;
  float *zeros;
  double mtheta1, mtheta2, ftheta1, ftheta2, lr1, lr2, lr3, lr4, cutlr2, cutlr3, cutlr4;
  double ld_pre=0, ld_post=0, le_pre=0, le_post=0, cutvol, denom, ppld;
  st_dim *mthetas, *fthetas;
 
  if ((zeros = malloc (sizeof (float) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((zero_idx = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);
  mthetas = &thetas->dims[0];
  fthetas = &thetas->dims[1];

  for (va = 0; va < dprimes->totalelems; va++) {

    /* See comment in calc_ppl_sexavg() */
    if (lr[va][0] == FLT_MAX)
      continue;

    if (va == zero_idx) {
      /* D' == 0 */
      for (vb = 1; vb < mthetas->numelems; vb++) {
	mtheta1 = mthetas->arr[vb - 1];
	mtheta2 = mthetas->arr[vb];
	for (vc = 1; vc < fthetas->numelems; vc++) {
	  ftheta1 = fthetas->arr[vc - 1];
	  ftheta2 = fthetas->arr[vc];
	  lr1 = lr[va][(vb - 1) * mthetas->numelems + vc - 1];
	  lr2 = lr[va][vb * mthetas->numelems + vc - 1];
	  lr3 = lr[va][(vb - 1) * mthetas->numelems + vc];
	  lr4 = lr[va][vb * mthetas->numelems + vc];
	  
	  if ((mtheta1 >= cutoff) || (ftheta1 >= cutoff)) {
	    /* entire region is outside cutoff area */
	    le_post += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	    
	  } else if ((mtheta2 <= cutoff) && (ftheta2 <= cutoff)) {
	    /* entire region inside cutoff area */
	    le_pre += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	    
	  } else if ((ftheta2 <= cutoff) && ((mtheta1 <= cutoff) && (cutoff <= mtheta2))) {
	    /* region straddles cutoff on male axis */
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    le_pre += (cutoff - mtheta1) * (ftheta2 - ftheta1) * (lr1 + cutlr2 + lr3 + cutlr4);
	    le_post += (mtheta2 - cutoff) * (ftheta2 - ftheta1) * (cutlr2 + lr2 + cutlr4 + lr4);
	    
	  } else if ((mtheta2 <= cutoff) && ((ftheta1 <= cutoff) && (cutoff <= ftheta2))) {
	    /* region straddles cutoff on female axis */
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutlr4 = lr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr4 - lr2);
	    le_pre += (cutoff - ftheta1) * (mtheta2 - mtheta1) * (lr1 + lr2 + cutlr3 + cutlr4);
	    le_post += (ftheta2 - cutoff) * (mtheta2 - mtheta1) * (cutlr3 + cutlr4 + lr3 + lr4);
	    
	  } else {
	    /* region straddles cutoff on both axes */
	    le_post += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    cutlr4 = cutlr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (cutlr4 - cutlr2);
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutvol = (cutoff - mtheta1) * (cutoff - ftheta1) * (lr1 + cutlr2 + cutlr3 + cutlr4);
	    le_post -= cutvol;
	    le_pre += cutvol;
	  }
	}
      }
    } else {
      /* D' != 0 */
      ndprimes++;
      for (vb = 1; vb < mthetas->numelems; vb++) {
	mtheta1 = mthetas->arr[vb - 1];
	mtheta2 = mthetas->arr[vb];
	for (vc = 1; vc < fthetas->numelems; vc++) {
	  ftheta1 = fthetas->arr[vc - 1];
	  ftheta2 = fthetas->arr[vc];
	  lr1 = lr[va][(vb - 1) * mthetas->numelems + vc - 1];
	  lr2 = lr[va][vb * mthetas->numelems + vc - 1];
	  lr3 = lr[va][(vb - 1) * mthetas->numelems + vc];
	  lr4 = lr[va][vb * mthetas->numelems + vc];
	  
	  if ((mtheta1 >= cutoff) || (ftheta1 >= cutoff)) {
	    /* entire region is outside cutoff area */
	    /*ld_post += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);*/
	    
	  } else if ((mtheta2 <= cutoff) && (ftheta2 <= cutoff)) {
	    /* entire region inside cutoff area */
	    ld_pre += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	    
	  } else if ((ftheta2 <= cutoff) && ((mtheta1 <= cutoff) && (cutoff <= mtheta2))) {
	    /* region straddles cutoff on male axis */
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    ld_pre += (cutoff - mtheta1) * (ftheta2 - ftheta1) * (lr1 + cutlr2 + lr3 + cutlr4);
	    /*ld_post += (mtheta2 - cutoff) * (ftheta2 - ftheta1) * (cutlr2 + lr2 + cutlr4 + lr4);*/
	    
	  } else if ((mtheta2 <= cutoff) && ((ftheta1 <= cutoff) && (cutoff <= ftheta2))) {
	    /* region straddles cutoff on female axis */
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutlr4 = lr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr4 - lr2);
	    ld_pre += (cutoff - ftheta1) * (mtheta2 - mtheta1) * (lr1 + lr2 + cutlr3 + cutlr4);
	    /*ld_post += (ftheta2 - cutoff) * (mtheta2 - mtheta1) * (cutlr3 + cutlr4 + lr3 + lr4);*/
	    
	  } else {
	    /* region straddles cutoff on both axes */
	    /*ld_post += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);*/
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    cutlr4 = cutlr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (cutlr4 - cutlr2);
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutvol = (cutoff - mtheta1) * (cutoff - ftheta1) * (lr1 + cutlr2 + cutlr3 + cutlr4);
	    /*ld_post -= cutvol;*/
	    ld_pre += cutvol;
	  }
	}
      }
    }
  }

  /* The 0.25 here is factored out of the calculation of area */
  /* We may be mis-using prior here by making no distiction between linkage prior and ld prior */
  ld_pre *= prior * (weight / (cutoff * cutoff)) * 0.25;
  ld_pre /= ndprimes;
  /*ld_post *= ((1 - weight) / (0.5 - cutoff)) * 0.25;*/
  /*ld_post /= ndprimes;*/
  le_pre *= (1 - prior) * (weight / (cutoff * cutoff)) * 0.25;
  le_post *= ((1 - weight) / (0.25 - (cutoff * cutoff))) * 0.25;
  if ((denom = ld_pre + ld_post + le_pre + le_post) == 0.0)
    return (0.0);
  ppld = (ld_pre + ld_post) / denom;
  return (ppld);
}


int parse_command_line (int argc, char **argv)
{
  int arg;

  opterr = 0;
  while ((arg = getopt (argc, argv, "sp:w:c:")) != -1) {
    if (arg == 's') {
      sexspecific = 1;
    } else if (arg == 'p') {
      prior = strtod (optarg, NULL);
    } else if (arg == 'w') {
      weight = strtod (optarg, NULL);
    } else if (arg == 'c') {
      cutoff = strtod (optarg, NULL);
    } else {
      fprintf (stderr, "bad option\n");
      exit (-1);
    }
  }
  return (optind);
}


int get_marker_line (st_marker *marker, FILE *fp)
{
  char buff[256];
  char *pa, *pb;

  if (fgets (buff, 256, fp) == NULL) {
    if (feof (fp))
      return (0);
    fprintf (stderr, "read error at line %d, %s\n", lineno, strerror (errno));
    exit (-1);
  }
  lineno++;

  if (((pa = strtok_r (buff, " \t\n", &pb)) == NULL) ||
      (strcmp (pa, "#") != 0))
    return (-1);
  
  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL)
    return (-1);
  if (((marker->num = (int) strtol (pa, NULL, 10)) == 0) && (errno != 0))
    return (-1);

  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL)
    return (-1);
  strcpy (marker->name1, pa);

  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL)
    return (-1);
  strcpy (marker->name2, pa);

  return (1);
}


int get_header_line (st_marker *marker, st_data *data, FILE *fp)
{ 
  char buff[256], token[32];
  char *pa, *pb, *pc=NULL;
  int actualcols, numlrcols=0;

  if (fgets (buff, 256, fp) == NULL) {
    if (feof (fp))
      return (0);
    fprintf (stderr, "read error at line %d, %s\n", lineno, strerror (errno));
    exit (-1);
  }
  lineno++;

  /* We'll only parse one header line per file, and require that the
   * format be consistant within the file. If this pointer is non-null,
   * we've parsed a header previsously, so just return.
   */
  if (marker->datacols != NULL)
    return (1);

  pa = strtok_r (buff, " \t\n", &pb);
  while (pa != NULL) {
    /*printf ("initial token '%s'\n", pa);*/
    strcpy (token, pa);
    while ((strchr (token, '(') != NULL) && (strrchr (token, ')') == NULL)) {
      if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL)
	return (-1);
      /*printf ("  tacking on '%s'\n", pa);*/
      strcat (token, pa);
    }
    /*printf ("final token '%s'\n", token);*/
    if (token[0] == '(')
      return (-1);

    if ((pa = strchr (token, '(')) == NULL) {
      actualcols = 1;
    } else {
      /* Some parenthesized expression */
      actualcols = 0;
      *pa = '\0';
      pa = strtok_r (++pa, ",", &pc);
      while (pa != NULL) {
	actualcols++;
	pa = strtok_r (NULL, ",", &pc);
      }
    }
    /*printf ("token is '%s', actualcols %d\n", token, actualcols);*/

    marker->numcols += actualcols;
    /*printf ("reallocating marker cols to %d\n", marker->numcols);*/
    if ((marker->datacols = realloc (marker->datacols, sizeof (int) * marker->numcols)) == NULL) {
      fprintf (stderr, "malloc failed, %s\n", strerror (errno));
      exit (-1);
    }

    if ((token[0] == 'D') && ((token[1] >= '0') && (token[1] <= '9')) &&
	((token[1] >= '0') && (token[1] <= '9'))) {
      /* A D-prime column */
      marker->numdprimes += actualcols;
      /*printf ("number of dprime cols %d\n", marker->numdprimes);*/
      for (; actualcols > 0; actualcols--) 
	marker->datacols[marker->numcols - actualcols] = DPRIME_COL;
      
    } else if (strcasecmp (token, "Theta") == 0) {
      /* A Theta column */
      if (sexspecific) {
	marker->numthetas += actualcols;
	/*printf ("number of theta cols %d\n", marker->numthetas);*/
	for (; actualcols > 0; actualcols--) 
	  marker->datacols[marker->numcols - actualcols] = THETA_COL;
      } else {
	marker->numthetas++;
	/*printf ("number of theta cols %d\n", marker->numthetas);*/
	marker->datacols[marker->numcols - actualcols] = THETA_COL;
	for (--actualcols; actualcols > 0; actualcols--) 
	  marker->datacols[marker->numcols - actualcols] = 0;
      }
      
    } else if ((strcasecmp (token, "AVG_LR") == 0) || (strcasecmp (token, "AVGLR") == 0)) {
      /* The average LR column */
      marker->datacols[marker->numcols - 1] = LR_COL;
      numlrcols++;
      /*printf ("number of lr cols %d\n", numlrcols);*/
      
    } else {
      /* Something else */
      marker->datacols[marker->numcols - 1] = 0;
    }

    /*printf ("\n");*/
    pa = strtok_r (NULL, " \t\n", &pb);
  }

  if (numlrcols != 1)
    return (-1);

  /* If there's no D' columns, then we've got a non-LD run, so we set that flag
   * We need at least one dimension of D's for our 2-dimensional array of avgLRs,
   * though, so we'll fake D' == 0 for every data line, and we need allocate
   * space for our dummy D'.
   */
  if (marker->numdprimes == 0) {
    marker->no_ld = 1;
    marker->numdprimes = 1;
  } else if (marker->numdprimes > 1) {
    marker->holey_grid = 1;
  }

  if (((data->dprimes = malloc (sizeof (float) * marker->numdprimes)) == NULL) ||
      ((data->thetas = malloc (sizeof (float) * marker->numthetas)) == NULL)) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  return (1);
}


int get_data_line (st_marker *marker, st_data *data, FILE *fp)
{
  char buff[256];
  char *pa, *pb;
  long start;
  int va = 0, dprimecnt=0, thetacnt=0;

  if ((start = ftell (fp)) == -1) {
    fprintf (stderr, "ftell failed, %s\n", strerror (errno));
    exit (-1);
  }
  if (fgets (buff, 256, fp) == NULL) {
    if (feof (fp))
      return (0);
    fprintf (stderr, "read error at line %d, %s\n", lineno, strerror (errno));
    exit (-1);
  }
  lineno++;
  
  if ((pa = strtok_r (buff, " (),\t\n", &pb)) == NULL) {
    fprintf (stderr, "unexpectedly short line %d\n", lineno);
    exit (-1);
  }

  /* If it's a new marker line, back the file pointer up, and return a 
   * non-error, non-dataline value.
   */
  if (strcmp (pa, "#") == 0) {
    if (fseek (fp, start, SEEK_SET) == -1) {
      fprintf (stderr, "fseek failed, %s\n", strerror (errno));
      exit (-1);
    }
    return (2);
  }

  while (1) {
    if (marker->datacols[va] == DPRIME_COL) {
      data->dprimes[dprimecnt++] = (float) strtod (pa, NULL);
    } else if (marker->datacols[va] == THETA_COL) {
      data->thetas[thetacnt++] = (float) strtod (pa, NULL);
    } else if (marker->datacols[va] == LR_COL) {
      data->lr = strtod (pa, NULL);
    }
    if (++va >= marker->numcols)
      break;
    if ((pa = strtok_r (NULL, " (),\t\n", &pb)) == NULL) {
      fprintf (stderr, "unexpectedly short line %d\n", lineno);
      exit (-1);
    }
  }

  /* If we've got a non-LD input file, we need to dummy up a D' of 0
   * for each data line, just so we've got a D' dimension for our
   * 2-dimensional array of avgLRs.
   */
  if (marker->no_ld) 
    data->dprimes[0] = 0;
  return (1);
}


int multi_insert (st_multidim *md, float *vals, int num)
{
  int va, idx, offset=0;
  st_dim *tmp;

  if (md->numdims != num) {
    if ((tmp = realloc (md->dims, sizeof (st_dim) * num)) == NULL)
      return (-1);
    md->dims = tmp;
    for (va = md->numdims; va < num; va++)
      memset (&(md->dims[va]), 0, sizeof (st_dim));
    md->numdims = num;
  }

  for (va = md->numdims - 1; va >= 0; va--) {
    if ((idx = insert (&(md->dims[va]), vals[va])) == -1)
      return (-1);
    if (va == md->numdims - 1)
      md->dims[va].dimsize = 1;
    else
      md->dims[va].dimsize = md->dims[va+1].numelems * md->dims[va+1].dimsize;
    offset += md->dims[va].dimsize * idx;
  }
  md->totalelems = md->dims[0].dimsize * md->dims[0].numelems;
  return (offset);
}


int multi_find (st_multidim *md, float *vals, int num)
{
  int va, idx, found, offset=0;

  if (md->numdims != num)
    return (-1);
  for (va = md->numdims - 1; va >= 0; va--) {
    idx = find (&(md->dims[va]), vals[va], &found);
    if (! found)
      return (-1);
    offset += idx * md->dims[va].dimsize;
  }
  return (offset);
}


int insert (st_dim *dim, float val)
{
  float *tmp;
  int idx, found;

  idx = find (dim, val, &found);
  if (found)
    return (idx);
  
  if (dim->arrsize < dim->numelems + 1) {
    if ((tmp = realloc (dim->arr, sizeof (float) * (dim->arrsize + 10))) == NULL)
      return (-1);
    dim->arr = tmp;
    dim->arrsize += 10;
  }
  if (idx < dim->numelems)
    memmove (dim->arr+idx+1, dim->arr+idx, sizeof (float) * (dim->numelems - idx));
  dim->arr[idx] = val;
  dim->numelems++;
  return (idx);
}


int find (st_dim *dim, float val, int *found)
{
  int va;

  *found = 0;
  if (dim->numelems == 0)
    return (dim->lastidx = 0);
  if (val == dim->arr[dim->lastidx]) {
    *found = 1;
    return (dim->lastidx);
  }

  if (val > dim->arr[dim->lastidx]) {
    for (va = dim->lastidx + 1; va < dim->numelems; va++) {
      if (val < dim->arr[va]) {
	return (dim->lastidx = va);
      } else if (val == dim->arr[va]) {
	*found = 1;
	return (dim->lastidx = va);
      }
    }
    return (dim->lastidx = dim->numelems);
    
  } else {   /*  val < dim->arr[dim->lastidx]  */
    for (va = 0; va <= dim->lastidx; va++) {
      if (val < dim->arr[va]) {
	return (dim->lastidx = va);
      } else if (val == dim->arr[va]) {
	*found = 1;
	return (dim->lastidx = va);
      }
    }
  }
  /* Should be impossible to reach this point */
  return (-1);
}


void multi_dump (st_multidim *md)
{
  int va, vb;

  for (va = 0; va < md->numdims; va++) {
    printf ("dim %d:", va);
    for (vb = 0; vb < md->dims[va].arrsize; vb++) {
      if (vb < md->dims[va].numelems) {
	printf (" %5.2f", md->dims[va].arr[vb]);
      } else {
	printf (" .....");
      }
    }
    printf (" (lastidx %d, dimsize %d)\n", md->dims[va].lastidx, md->dims[va].dimsize);
  }
  printf ("totalelems %d\n", md->totalelems);
  return;
}
