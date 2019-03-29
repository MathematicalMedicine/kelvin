/* Random Access DAta Storage for MLip Mutant Version (RADSMM) */

/* these routines Do Input/Output on an already opened file. */

#include <stddef.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <values.h>
#include "RADSMM.h"

/* This is an internal routine that picks the right fp and
moves it the position for a read or write. Returns the result of
lseek */

int RADSMM_seek( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
   long genefreq_index, long penet_index, long qmodel_index, long diseq_index, int *fpndx )


{
	long long 	offset;  /* non-standard */
	off_t		stat, short_offset;


        switch (header->model_type) {
            case 'D':
                qmodel_index = 0;
                break;
            case 'Q':
                penet_index = 0;
                genefreq_index = 0;
                break;
        }
 
        switch (header->marker_type) {
            case 'M':
                theta_index = 0;
        }
        switch (header->use_Diseq) {
            case 'N':
                diseq_index = 0;
        }


	if ( header->debug_level & 28 ) {  /* if debug flag set, do range checks */
	   if ( (ped_index >= header->pedigree_count) || (ped_index < 0)  ) {
		printf( " RADSMM:Invalid seek pedigree index %li\n",ped_index);
		return RADSMM_ERROR_badparam;
	   }
	   stat = header->theta_count;
	   if ( header->theta_matrix_type == 'G' ) stat = stat * header->theta_count;
	   if ( (theta_index >= stat) || (theta_index < 0) ) {
		printf( " RADSMM:Invalid seek theta index %li\n",theta_index);
		return RADSMM_ERROR_badparam;
	   }
	   if ( (genefreq_index >= header->geneFreq_count) || (genefreq_index < 0) ) {
		printf( " RADSMM:Invalid seek gene freq index %li\n",genefreq_index);
		return RADSMM_ERROR_badparam;
	   }
	   if ( (marker_index >= header->marker_count) || (marker_index < 0) ) {
		printf( " RADSMM:Invalid seek marker index %li\n",marker_index);
		return RADSMM_ERROR_badparam;
	   }
           if ( (penet_index >= header->penetrance_count) || (penet_index < 0) ) {
	        printf( " RADSMM:Invalid seek penetrance index %li\n",penet_index);
	        return RADSMM_ERROR_badparam;
	   }

	   if ( (diseq_index >= header->diseq_count) || (diseq_index < 0) ) {
	        printf( " RADSMM:Invalid seek diseq index %li %li\n",diseq_index,header->diseq_count);
	        return RADSMM_ERROR_badparam;
	   }

	   if ( (qmodel_index >= header->qmodel_count) || (qmodel_index < 0) ) {
	        printf( " RADSMM:Invalid seek qmodel index %li\n",qmodel_index);
	        return RADSMM_ERROR_badparam;
	   }

	}

/*  Taken straight from 'doc'
 A: Pedigree, Marker, MaleTheta, FemaleTheta, GeneFreq, Penetrance, Qmodel, diseq
 B: Pedigree, GeneFreq, Penetrance, Qmodel, MaleTheta, FemaleTheta, diseq, Marker
 C: Marker, Pedigree, Qmodel, MaleTheta, FemaleTheta, GeneFreq, Penetrance, diseq
 D: Pedigree, MaleTheta, FemaleTheta, Marker, Genefreq, Penetrance, qmodel, diseq
 E: Pedigree, qmodel, GeneFreq, Penetrance, diseq, Marker, MaleTheta, FemaleTheta
 F: Marker, GeneFreq, Penetrance, Qmodel, MaleTheta, FemaleTheta, diseq, Pedigree
*/


	switch (header->ordering) {
	   case 'A':
	      offset = ped_index;
	      offset = offset*header->marker_count + marker_index;
	      if ( header->theta_matrix_type == 'D' ) {
	         offset = offset*header->theta_count + theta_index;
	      } else {  /* just assume it's G */
 	         offset = offset*(header->theta_count*header->theta_count) + theta_index;
	      }
	      offset = offset*header->geneFreq_count + genefreq_index;
	      offset = offset*header->penetrance_count + penet_index;
              offset = offset*header->qmodel_count + qmodel_index;
              offset = offset*header->diseq_count + diseq_index;
	      break;
	   case 'B':
	      offset = ped_index;
	      offset = offset*header->geneFreq_count + genefreq_index;
	      offset = offset*header->penetrance_count + penet_index;
              offset = offset*header->qmodel_count + qmodel_index;
	      if ( header->theta_matrix_type == 'D' ) {
	         offset = offset*header->theta_count + theta_index;
	      } else {  /* just assume it's G */
	         offset = offset*(header->theta_count*header->theta_count) + theta_index;
	      }
              offset = offset*header->diseq_count + diseq_index;
	      offset = offset*header->marker_count + marker_index;
	      break;
	   case 'C':
	      offset = marker_index;
	      offset = offset*header->pedigree_count + ped_index;
              offset = offset*header->qmodel_count + qmodel_index;
	      if ( header->theta_matrix_type == 'D' ) {
	         offset = offset*header->theta_count + theta_index;
	      } else {  /* just assume it's G */
 	         offset = offset*(header->theta_count*header->theta_count) + theta_index;
	      }
	      offset = offset*header->geneFreq_count + genefreq_index;
	      offset = offset*header->penetrance_count + penet_index;
              offset = offset*header->diseq_count + diseq_index;
	      break;
	   case 'D':
	      offset = ped_index;
	      if ( header->theta_matrix_type == 'D' ) {
	         offset = offset*header->theta_count + theta_index;
	      } else {  /* just assume it's G */
 	         offset = offset*(header->theta_count*header->theta_count) + theta_index;
	      }
	      offset = offset*header->marker_count + marker_index;
	      offset = offset*header->geneFreq_count + genefreq_index;
	      offset = offset*header->penetrance_count + penet_index;
              offset = offset*header->qmodel_count + qmodel_index;
              offset = offset*header->diseq_count + diseq_index;
	      break;
	   case 'E':
	      offset = ped_index;
              offset = offset*header->qmodel_count + qmodel_index;
	      offset = offset*header->geneFreq_count + genefreq_index;
	      offset = offset*header->penetrance_count + penet_index;
              offset = offset*header->diseq_count + diseq_index;
	      offset = offset*header->marker_count + marker_index;
	      if ( header->theta_matrix_type == 'D' ) {
	         offset = offset*header->theta_count + theta_index;
	      } else {  /* just assume it's G */
 	         offset = offset*(header->theta_count*header->theta_count) + theta_index;
	      }
	      break;
	   case 'F':
	      offset = marker_index;
	      offset = offset*header->geneFreq_count + genefreq_index;
	      offset = offset*header->penetrance_count + penet_index;
              offset = offset*header->qmodel_count + qmodel_index;
	      offset = offset*header->marker_count + marker_index;
	      if ( header->theta_matrix_type == 'D' ) {
	         offset = offset*header->theta_count + theta_index;
	      } else {  /* just assume it's G */
 	         offset = offset*(header->theta_count*header->theta_count) + theta_index;
	      }
              offset = offset*header->diseq_count + diseq_index;
	      offset = offset*header->pedigree_count + ped_index;
	      break;
	   default:
	      if ( header->debug_level & 28 ) {
	         printf( " RADSMM:Error BAD Ordering	'%c'\n",header->ordering);
	      }
	      return RADSMM_ERROR_badparam;
	}

	if ( offset < 0 ) {
	  if( header->debug_level & 63 ) 
	    printf(" RADSMM:seeking before start of file\n");
	  return RADSMM_ERROR_badparam;
	}

	if ( header->number_of_files >= 1 ) {
		/* more than 1 file */
		*fpndx = 1 + ( offset / header->chunks_per_file );
		short_offset = offset - ((*fpndx-1)* header->chunks_per_file);
		short_offset = short_offset * header->chunk_size;
/*		if ( header->debug_level & 16 ) {
		   printf( " offset=%lli, fpndx=%i, short=%li",offset,*fpndx,short_offset);
		   printf( " chunks per file =%li , chunk_size=%i\n",header->chunks_per_file,header->chunk_size);
		} */
	} else {
		/* only 1 file */
	    *fpndx = 0;
	    offset = offset * header->chunk_size;
	    short_offset = offset + header->start_of_data;
	}
 
	stat = lseek( header->fp[*fpndx], short_offset, SEEK_SET );
	if ( (stat == -1 ) && ( (header->debug_level & 3 ) != 0 ) ){
	    printf( " RADSMM:* seek error: err#:%i, offset=%lli,%li,%i\n",errno,offset,short_offset,*fpndx);
	    return RADSMM_ERROR_lseek;
	}
/*	if ( header->debug_level & 16 ) {
		printf("*** seeked to %i: %i+ %i:%i:%i:%i:%i\n", stat, 
		 header->start_of_data,ped_index,theta_index,penet_index,
		 genefreq_index,marker_index); 
	} */

	return 0;
}

                 
