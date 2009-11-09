
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef enum
{ DAD = 0, MOM = 1 } PARENTS;
typedef enum
{ MALE = 1, FEMALE = 2 } SEX;

/* maximum length for pedigree label, individual IDs */
#define MAX_PED_LABEL_LEN       128

#define ND NORMAL_DISTRIBUTION	/* normal distribution */
#define TD T_DISTRIBUTION	/* t distribution */

#include "pedigree.h"
#include "locus.h"
#include "../config/config.h"

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
void adjustQuantitativeTraits (PedigreeSet *pPedigreeSet);
void getPedigreeTraitRange (PedigreeSet *pPedigreeSet, double *min, double *max);
void getPedigreeSampleStdev (PedigreeSet *pPedigreeSet, double *mean, double *stdev);
void renumberLiabilityClasses (PedigreeSet *pPedigreeSet);

extern char *flexBuffer;
extern int flexBufferSize;

#endif
