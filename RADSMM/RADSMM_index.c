/* RADSMM */

/* This are index routines, they return the location of
a value in list of values. */

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "RADSMM.h"

long RADSMM_index_pedigree( RADSMM_header_type *header, int pedid )
{
	long	i;
	RADSMM_pedigree_type *pedptr;

	
	pedptr = header->pedigree_list;

	for (i=0; i<header->pedigree_count; i++ )
	{
		if ( *pedptr == pedid ) return i;
		pedptr++;
	}
	return RADSMM_ERROR_value_not_in_list;
}

long RADSMM_index_theta( RADSMM_header_type *header, double male_theta, 
						     double female_theta )
{
	long	i, mindex, findex;
	double	close, closev;
	RADSMM_theta_type *tptr;

	if ( header->model_type == 'Q' ) { 
		return RADSMM_ERROR_wrong_model;
	}

	if ( header->theta_matrix_type == 'D' )
	{
		if ( fabs(male_theta-female_theta) > 2.0E-5 ) {
			if (header->debug_level & 3 )
			    printf(" Male and Female thetas not equal on Theta Type = D %f %f\n",male_theta,female_theta);
			return RADSMM_ERROR_badparam;
		}
		tptr = header->theta_list;
		for (i=0; i < header->theta_count; i++ )
		{
			if ( fabs(*tptr-male_theta) < 1.0E-5 ) return i;
			tptr++;
		}
		if ( header->debug_level & 3 ) {
		    printf(" Theta value not found in list %10.8f\n",male_theta);
		    closev = close = 9e99;
		    tptr = header->theta_list;
		    for (i=0; i < header->theta_count; i++ ) {
			if ( fabs(*tptr-male_theta) < close ) {
				close = fabs(*tptr-male_theta);
				closev = *tptr;
			}
			tptr++;
		    }
		    printf(" closest value found is %10.8f with an offset of %e",
				closev, close );
		}
		return RADSMM_ERROR_value_not_in_list;

	}
	else if ( header->theta_matrix_type == 'G' )
	{
		mindex = findex = -1;
		tptr = header->theta_list;
		for (i=0; i < header->theta_count; i++ )
		{
			if ( fabs(*tptr-male_theta) < 1.0E-5 ) {
				mindex = i;
				break;
			}
			tptr++;
		}
		if ( mindex == -1 ) {
			if ( header->debug_level & 3 )
			    printf(" Male Theta value not found in list %10.8f\n",male_theta);
			return RADSMM_ERROR_value_not_in_list;
		}
		tptr = header->theta_list;
		for (i=0; i < header->theta_count; i++ )
		{
			if ( fabs(*tptr-female_theta) < 1.0E-5 ) {
				findex = i;
				break;
			}
			tptr++;
		}
		if ( findex == -1 )  {
			if ( header->debug_level & 3 )
			    printf(" Female Theta value not found in list %10.8f\n",female_theta);
			return RADSMM_ERROR_value_not_in_list;
		}
		
		return mindex*header->theta_count + findex;
	}
	return RADSMM_ERROR_internal;	
}

long RADSMM_index_penetrance( RADSMM_header_type *header,  int LC,
			      float pen1, float pen2, float pen3  )
{
	long	i;
	RADSMM_penetrance_type *penptr;


	if ( header->model_type == 'Q' ) { 
		return RADSMM_ERROR_wrong_model;
	}

        if ( (LC < 0) || (LC >= header->LC_count) ) {
	   if ( header->debug_level & 3 ) {
		printf(" ** RADSMM_index_penet: Bad LC value :%i\n",LC ); 
	   }
	   return RADSMM_ERROR_badparam;
        }

	penptr = header->penetrance_list[LC];

	for (i=0; i<header->penetrance_count; i++ )
	{
		if ( ( fabsf(penptr->penetrance[0]-pen1) < 1.0E-5 ) &&
		     ( fabsf(penptr->penetrance[1]-pen2) < 1.0E-5 ) &&
		     ( fabsf(penptr->penetrance[2]-pen3) < 1.0E-5 ) ) return i;
		penptr++;
	}
	if ( header->debug_level & 3 ) {
		printf(" ** no such pen combination:%f,%f,%f\n",pen1,pen2,pen3 ); 
	}
	return RADSMM_ERROR_value_not_in_list;
}

long RADSMM_index_genefreq( RADSMM_header_type *header, double gf )
{
	long	i;
	RADSMM_geneFreq_type *gfptr;

	if ( header->model_type == 'Q' ) { 
		return RADSMM_ERROR_wrong_model;
	}

	gfptr = header->geneFreq_list;

	for (i=0; i<header->geneFreq_count; i++ )
	{
		if ( fabs(*gfptr-gf) < 1.0E-6 ) return i;
		gfptr++;
	}
	if ( header->debug_level & 24 ) {
		printf(" ** no such genefreq:%f\n",gf ); 
	}
	return RADSMM_ERROR_value_not_in_list;
}


/* long RADSMM_index_qmodel( RADSMM_header_type *header,
			      float m1,
			      float m2,
			      float m3  )
{
	long	i;
	mean_type *ptr;


	if ( header->model_type != 'Q' ) { 
		return RADSMM_ERROR_wrong_model;
	}
	ptr = header->mean_list;

	for (i=0; i<header->mean_count; i++ )
	{
		if ( ( fabsf(ptr->means[0]-m1) < 1.0E-5 ) &&
		     ( fabsf(ptr->means[1]-m2) < 1.0E-5 ) &&
		     ( fabsf(ptr->means[2]-m3) < 1.0E-5 ) ) return i;
		ptr++;
	}
	if ( header->debug_level & 3 ) {
		printf(" ** no such means combination:%f,%f,%f\n",m1,m2,m3 ); 
	}
	return RADSMM_ERROR_value_not_in_list;
}
*/

long RADSMM_index_diseq( RADSMM_header_type *header, double lam )
{
	long	i;
	RADSMM_diseq_type *ptr;


	if ( header->use_Diseq == 'N' ) { 
		return RADSMM_ERROR_wrong_model;
	}

	ptr = header->diseq_list;

	for (i=0; i<header->diseq_count; i++ )
	{
		if ( fabs(*ptr-lam) < 1.0E-6 ) return i;
		ptr++;
	}
	if ( header->debug_level & 24 ) {
		printf(" ** no such diseq value:%f\n",lam ); 
	}
	return RADSMM_ERROR_value_not_in_list;
}


