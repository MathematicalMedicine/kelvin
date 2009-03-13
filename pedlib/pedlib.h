
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __PEDLIB_H__
#define __PEDLIB_H__

/********************************************************************
 * Include file for Kelvin pedigree components.
 * This library reads in pedigree file and build pedigree structure
 * reads in marker map file and build the initial marker map
 * does genotype elimination, allele set recoding
 * provides functionality and API for likelihood calculation
 *******************************************************************/

/* Analysis mode (LE or LD) */
#define LINKAGE_EQUILIBRIUM		0
#define LINKAGE_DISEQUILIBRIUM	1

/* Analysis type in a different perspective: 2point or multipoint */
#define TWOPOINT	0
#define MULTIPOINT	1

/* Trait value type */
#define DICHOTOMOUS	0

/* The main difference between DT and QT are in the way how conditional
 * penetrances are determined. For DT, penetrances are given, while for
 * QT, they are calculated using means and variances. */
#define QUANTITATIVE	1

/* Combination of dichotomous and quantitative trait 
 * some individuals in the pedigree has known quantitative trait values
 * individuals with unknown trait values are handled in a way similar to 
 * under dichotomouse trait (using some cutoff values to calculate conditional
 * penetrance given genotypes with CDF) */
#define COMBINED	2

/* Marker to marker analysis types. FALSE means no analysis. */
#define MARKERTOMARKER 1
#define ADJACENTMARKER 2

typedef enum
{ DAD = 0, MOM = 1 } PARENTS;
typedef enum
{ MALE = 1, FEMALE = 2 } SEX;

/* maximum length for pedigree label, individual IDs */
#define MAX_PED_LABEL_LEN       128

#define TP TWOPOINT		/* 2 point */
#define MP MULTIPOINT		/* multipoint */
#define SA SEX_AVERAGED
#define SS SEX_SPECIFIC
#define DT DICHOTOMOUS		/* dichotomous trait */
#define QT QUANTITATIVE		/* quantitative trait */
#define CT COMBINED		/* combined dichotomous/quantitative trait */
#define MM MARKERTOMARKER	/* marker to marker analysis type */
#define AM ADJACENTMARKER
#define ND NORMAL_DISTRIBUTION	/* normal distribution */
#define TD T_DISTRIBUTION	/* t distribution */
#define NPENET(x) ((x)*(x))

#include "pedigree.h"
#include "locus.h"
#include "model_options.h"
extern struct ModelOptions modelOptions;
#include "model_type.h"
extern struct ModelType modelType;
#include "model_range.h"
extern struct ModelRange modelRange;

/* pedigree related function prototypes */
int read_pedfile (char *sPedfileName, PedigreeSet * pedigreeSet);
int read_ccfile (char *ccFileName, PedigreeSet * pedigreeSet);
Person *find_person (Pedigree * pPed, char *sPersonID);
Pedigree *find_pedigree (PedigreeSet * pPedSet, char *sPedID);

int read_mapfile (char *sMapfileName);
int read_datafile (char *sMarkerfileName);

int allele_set_recoding (int locus, Pedigree * pPedigree);
int pedigree_genotype_elimination (int locus, Pedigree * pPedigree);
int compute_likelihood (PedigreeSet * pPedigreeList);
void pedigreeSetPolynomialClearance (PedigreeSet * pPedigreeList);
int compute_pedigree_likelihood (Pedigree * pPedigree);
void *MALLOC (char *description, size_t size);
void FREE (char *description, void *ptr);
void *REALLOC (char *description, void *ptr, size_t size);

extern char *flexBuffer;
extern int flexBufferSize;

#endif
