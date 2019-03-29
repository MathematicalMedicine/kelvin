/* Random Access DAta Storage for MLip  (RADSMM) */


/* these routines Do Input/Output on an already opened file. */

#include <stddef.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <values.h>
#include "RADSMM.h"

int RADSMM_read_data( RADSMM_header_type *header, 
                        long ped_index,
                        long marker_index,
                        long theta_index,
                        long genefreq_index,
                        long penet_index,
                        long qmodel_index, 
                        long diseq_index,
			double *data )
{
   int	fpndx;
   long	rc;
   ssize_t	stat;
   double  dbuffer;
   float   fbuffer;

	if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, 
             genefreq_index, penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
	    return rc;
	}

/*	printf(" %i %s = %i\n", header->elements_per_item, 
		&header->element_data_type, header->chunk_size); */
	if ( header->element_data_type == 'D' ) {
		stat = read( header->fp[fpndx], &dbuffer, sizeof(double) );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readd %i, errno=%i, %e\n",stat,errno,dbuffer);
			return RADSMM_ERROR_reading;
		}
		*data = dbuffer;
/*		if ( header->debug_level & 16) 
		  printf( " *** readin out %E(%E) \n",dbuffer,*data ); */
	} else if ( header->element_data_type == 'F' ) {
		stat = read( header->fp[fpndx], &fbuffer, sizeof(float) );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readf %i, errno=%i, %E\n",stat,errno,fbuffer);
			return RADSMM_ERROR_reading;
		}
		*data = fbuffer;
/*		if ( (header->debug_level & 16) !=0) 
		  printf( " *** readin out %E(%E) \n",fbuffer,data ); */
	} else
		return RADSMM_ERROR_badparam;

	return 0;
}

int RADSMM_read_float( RADSMM_header_type *header, 
                        long ped_index,
                        long marker_index,
                        long theta_index,
                        long genefreq_index,
                        long penet_index,
                        long qmodel_index, 
                        long diseq_index,
			float *data )
{
	int	fpndx;
	long	rc;
	ssize_t	stat;
        double  dbuffer;
        float   fbuffer;

	if ( header->open_flag != 'O' ) return RADSMM_ERROR_not_open;

	if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, 
             genefreq_index, penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
		return rc;
	}

/*	printf(" %i %s = %i\n", header->elements_per_item, 
		&header->element_data_type, header->chunk_size); */
	if ( header->element_data_type == 'D' ) {
		stat = read( header->fp[fpndx], &dbuffer, sizeof(double) );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readd %i, errno=%i, %e\n",stat,errno,dbuffer);
			return RADSMM_ERROR_reading;
		}
                if ( (dbuffer > MAXFLOAT) || (dbuffer < MINFLOAT) ) {
                        return RADSMM_ERROR_outofrange;
                    }
		*data = dbuffer;
/*		if ( header->debug_level & 16) 
		  printf( " *** readin out %E(%E) \n",dbuffer,*data ); */
	} else if ( header->element_data_type == 'F' ) {
		stat = read( header->fp[fpndx], &fbuffer, sizeof(float) );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readf %i, errno=%i, %E\n",stat,errno,fbuffer);
			return RADSMM_ERROR_reading;
		}
		*data = fbuffer;
/*		if ( (header->debug_level & 16) !=0) 
		  printf( " *** readin out %E(%E) \n",fbuffer,data ); */
	} else
		return RADSMM_ERROR_badparam;

	return 0;
}

int RADSMM_read_list_float( RADSMM_header_type *header, 
                        long ped_index,
                        long marker_index,
                        long theta_index,
                        long genefreq_index,
                        long penet_index,
                        long qmodel_index, 
                        long diseq_index,
			float data[],
                        long ndxcnt )
{
	int	i, fpndx;
	long	rc;
	ssize_t	stat;
        double  dbuffer[ndxcnt];


	if ( header->open_flag != 'O' ) return RADSMM_ERROR_not_open;

        i = RADSMM_range_check(  header, ped_index, marker_index, theta_index, 
                 genefreq_index, penet_index, qmodel_index, diseq_index, ndxcnt );
        if ( i != RADSMM_ERROR_success ) {
           if ( header->debug_level & 28 )
                printf( " ndxcnt out of range %li\n",ndxcnt);
           return i;
        }

	if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, 
             genefreq_index, penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
		return rc;
	}



/*	printf(" %i %s = %i\n", header->elements_per_item, 
		&header->element_data_type, header->chunk_size); */
	if ( header->element_data_type == 'D' ) {
		stat = read( header->fp[fpndx], dbuffer, sizeof(double)*ndxcnt );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readd %i, errno=%i, %e\n",stat,errno,dbuffer[0]);
			return RADSMM_ERROR_reading;
		}
                for ( i=0 ; i<ndxcnt ; i++ ) {
                    if ( (data[i] > MAXFLOAT) || (data[i] < MINFLOAT) ) {
                        return RADSMM_ERROR_outofrange;
                    }
		   data[i] = dbuffer[i];
		}
/*		if ( header->debug_level & 16) 
		  printf( " *** readin out %E(%E) \n",dbuffer,*data ); */
	} else if ( header->element_data_type == 'F' ) {
		stat = read( header->fp[fpndx], data, sizeof(float)*ndxcnt );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readf %i, errno=%i, %E\n",stat,errno,data[0]);
			return RADSMM_ERROR_reading;
		}
/*		if ( (header->debug_level & 16) !=0) 
		  printf( " *** readin out %E(%E) \n",fbuffer,data ); */
	} else
		return RADSMM_ERROR_badparam;

	return 0;
}

int RADSMM_read_list_double( RADSMM_header_type *header, 
                        long ped_index,
                        long marker_index,
                        long theta_index,
                        long genefreq_index,
                        long penet_index,
                        long qmodel_index, 
                        long diseq_index,
			double data[],
                        long ndxcnt )
{
	int	i, fpndx;
	long	rc;
	ssize_t	stat;
        float  fbuffer[ndxcnt];


	if ( header->open_flag != 'O' ) return RADSMM_ERROR_not_open;

        i = RADSMM_range_check(  header, ped_index, marker_index, theta_index, 
                 genefreq_index, penet_index, qmodel_index, diseq_index, ndxcnt );
        if ( i != RADSMM_ERROR_success ) {
           if ( header->debug_level & 28 )
                printf( " ndxcnt out of range %li\n",ndxcnt);
           return i;
        }

	if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, 
             genefreq_index, penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
		return rc;
	}


/*	printf(" %i %s = %i\n", header->elements_per_item, 
		&header->element_data_type, header->chunk_size); */
	if ( header->element_data_type == 'D' ) {
		stat = read( header->fp[fpndx], data, sizeof(double)*ndxcnt );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readd %i, errno=%i, %e\n",stat,errno,data[0]);
			return RADSMM_ERROR_reading;
		}
/*		if ( header->debug_level & 16) 
		  printf( " *** readin out %E(%E) \n",dbuffer,*data ); */
	} else if ( header->element_data_type == 'F' ) {
		stat = read( header->fp[fpndx], fbuffer, sizeof(float)*ndxcnt );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err readf %i, errno=%i, %E\n",stat,errno,fbuffer[0]);
			return RADSMM_ERROR_reading;
		}
                for ( i=0 ; i<ndxcnt ; i++ )
		   data[i] = fbuffer[i];
/*		if ( (header->debug_level & 16) !=0) 
		  printf( " *** readin out %E(%E) \n",fbuffer,data ); */
	} else
		return RADSMM_ERROR_badparam;

	return 0;
}

