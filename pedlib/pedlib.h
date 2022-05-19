#ifndef __PEDLIB_H__
#define __PEDLIB_H__
/* Copyright (C) 2006, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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


#define ND NORMAL_DISTRIBUTION	/* normal distribution */
#define TD T_DISTRIBUTION	/* t distribution */

#include "pedigree.h"
#include "locus.h"
#include "../config/model.h"

/* pedigree related function prototypes */
int read_pedfile (char *sPedfileName, PedigreeSet * pedigreeSet);
int read_ccfile (char *ccFileName, PedigreeSet * pedigreeSet);
Person *find_person (Pedigree * pPed, char *sPersonID);
Pedigree *find_pedigree (PedigreeSet * pPedSet, char *sPedID);

int read_mapfile (char *sMapfileName, int sexspecific);
int read_datafile (char *sMarkerfileName);

int allele_set_recoding (int locus, Pedigree * pPedigree);
int pedigree_genotype_elimination (int locus, Pedigree * pPedigree);
void pedigreeSetPolynomialClearance (PedigreeSet * pPedigreeList);
int compute_pedigree_likelihood (Pedigree * pPedigree);
void adjustQuantitativeTraits (PedigreeSet *pPedigreeSet);
int checkQtTraitRanges (PedigreeSet *pPedigreeSet);
void getPedigreeSampleStdev (PedigreeSet *pPedigreeSet, double *mean, double *stdev);
void renumberLiabilityClasses (PedigreeSet *pPedigreeSet);

/* Genotype Elimination State - this is to help get rid off some warning messages 
 * initialize state (0) - first state - genotype elimination warnings should be on
 * running state (1) - second state - genotype elimination warnings should be off
 */
int genoElimState;

#endif
