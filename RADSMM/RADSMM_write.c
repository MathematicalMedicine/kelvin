/* Random Access DAta Storage for MLip  (RADSML) */


/* these routines Do Input/Output on an already opened file. */

#include <stddef.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <values.h>
#include "RADSMM.h"


int RADSMM_range_check(   RADSMM_header_type *header,
                        long ped_index,
                        long marker_index,
                        long theta_index,
                        long genefreq_index,
                        long penet_index,
                        long qmodel_index,
                        long diseq_index,
                        long ndxcnt )
{

	if ( ndxcnt <= 0 ) return RADSMM_ERROR_badparam;
        switch ( header->ordering ) {
           case 'A':
           case 'C':
           case 'D':
              if ( (diseq_index+ndxcnt) > header->diseq_count )
                return RADSMM_ERROR_outofrange;
              break;
           case 'B':
              if ( (marker_index+ndxcnt) > header->marker_count )
                return RADSMM_ERROR_outofrange;
              break;
           case 'E':
              if ( header->theta_matrix_type == 'D') {
                 if ( (theta_index+ndxcnt) > header->theta_count )
                    return RADSMM_ERROR_outofrange;
              } else {
                 if ( (theta_index+ndxcnt) > (header->theta_count*header->theta_count) )
                    return RADSMM_ERROR_outofrange;
              }  
              break;
        }

	return RADSMM_ERROR_success;
}

int RADSMM_write_data( RADSMM_header_type *header, 
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

	if ( header->open_flag != 'O' ) {
	    if ( header->debug_level & 28 ) printf( " *** File not open %c\n",header->open_flag);
	    return RADSMM_ERROR_not_open;
	}

	if ( header->debug_level & 128 ) {
/*  Lets do a read check */
	   if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, genefreq_index, 
                        penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
		return rc;
	   }
	   if ( header->element_data_type == 'D' ) {
		stat = read( header->fp[fpndx], &dbuffer, sizeof(double) );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err chk readd %i, errno=%i, %e\n",stat,errno,dbuffer);
			return RADSMM_ERROR_reading;
		}
		if ( (fabs(dbuffer - RADSMM_EMPTY ) < 1.0E33 ) &&
                     (fabs(dbuffer - RADSMM_IGNORED ) < 1.0E33 ) &&
                     (fabs(dbuffer - *data) < fabs(*data/5.0E5) ) ) 
		 	    return RADSMM_ERROR_writeover_valid_data;

	   } else if ( header->element_data_type == 'F' ) {
		stat = read( header->fp[fpndx], &fbuffer, sizeof(float) );
		if ( stat == -1 ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err chk readf %i, errno=%i, %E\n",stat,errno,fbuffer);
			return RADSMM_ERROR_reading;
		}
		if ( (fabs(fbuffer - RADSMM_EMPTY ) < 1.0E33 ) &&
                     (fabs(fbuffer - RADSMM_IGNORED ) < 1.0E33 ) &&
                     (fabs(fbuffer - *data) < fabs(*data/5.0E5) ) )
		       	         return RADSMM_ERROR_writeover_valid_data;
	   }
	} /* if debug check */

	if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, genefreq_index, 
                        penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
		return rc;
	}
	if ( header->element_data_type == 'D' ) {
		dbuffer = *data;
		stat = write( header->fp[fpndx], &dbuffer, sizeof( double ) );
		if ( stat != sizeof(double) ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err writed %i, errno=%i\n",stat,errno);
			return RADSMM_ERROR_writing;
		}
/*		if ( header->debug_level & 16 ) 
		  printf( " *** wrote out %E(%E) \n",dbuffer,*data);  */
	} else if ( header->element_data_type == 'F' ) {
		fbuffer = *data;
		stat = write( header->fp[fpndx], &fbuffer, sizeof( float ) );
		if ( stat != sizeof(float) ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err writef %i, errno=%i\n",stat,errno);
			return RADSMM_ERROR_writing;
		}
/*		if ( header->debug_level & 16 ) 
		  printf( " *** wrote out %E(%E) \n",fbuffer,*data); */
	} else
		return RADSMM_ERROR_badparam;

	return 0;
}

int RADSMM_write_list_float(  RADSMM_header_type *header, 
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


	i = RADSMM_range_check( header, ped_index, marker_index, theta_index, genefreq_index, 
                                 penet_index, qmodel_index, diseq_index, ndxcnt );
	if ( i != RADSMM_ERROR_success ) {
	   if ( header->debug_level & 28 )
		printf( " ndxcnt out of range %li\n",ndxcnt);
	   return i;
	}

	if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, genefreq_index, 
                        penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
		return rc;
	}
	if ( header->element_data_type == 'D' ) {
                for ( i=0 ; i< ndxcnt ; i++ )
		   dbuffer[i] = data[i];
		stat = write( header->fp[fpndx], dbuffer, sizeof(double)*ndxcnt );
		if ( stat != sizeof(double)*ndxcnt ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err writed %i, errno=%i\n",stat,errno);
			return RADSMM_ERROR_writing;
		}
/*		if ( header->debug_level & 16 ) 
		  printf( " *** wrote out %E(%E) \n",dbuffer,*data);  */
	} else if ( header->element_data_type == 'F' ) {
		stat = write( header->fp[fpndx], data, sizeof(float)*ndxcnt );
		if ( stat != sizeof(float)*ndxcnt ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err writef %i, errno=%i\n",stat,errno);
			return RADSMM_ERROR_writing;
		}
/*		if ( header->debug_level & 16 ) 
		  printf( " *** wrote out %E(%E) \n",fbuffer,*data); */
	} else
		return RADSMM_ERROR_badparam;

	return 0;
}

int RADSMM_write_list_double(  RADSMM_header_type *header, 
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
        double  fbuffer[ndxcnt];

	if ( header->open_flag != 'O' ) return RADSMM_ERROR_not_open;

	i = RADSMM_range_check( header, ped_index, marker_index, theta_index, genefreq_index, 
                                 penet_index, qmodel_index, diseq_index, ndxcnt );
	if ( i != RADSMM_ERROR_success ) {
	   if ( header->debug_level & 28 )
		printf( " ndxcnt out of range %li\n",ndxcnt);
	   return i;
	}


	if ( (rc=RADSMM_seek( header, ped_index, marker_index, theta_index, genefreq_index, 
                        penet_index, qmodel_index, diseq_index, &fpndx )) < 0 ) {
		return rc;
	}
	if ( header->element_data_type == 'D' ) {
		stat = write( header->fp[fpndx], data, sizeof(double)*ndxcnt );
		if ( stat != sizeof(double)*ndxcnt ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err writed %i, errno=%i\n",stat,errno);
			return RADSMM_ERROR_writing;
		}
/*		if ( header->debug_level & 16 ) 
		  printf( " *** wrote out %E(%E) \n",dbuffer,*data);  */
	} else if ( header->element_data_type == 'F' ) {
                for ( i=0 ; i< ndxcnt ; i++ ) {
		    if ( (data[i] > MAXFLOAT) || (data[i] < MINFLOAT) ) {
			return RADSMM_ERROR_outofrange;
		    } 
		    fbuffer[i] = data[i];
                }
		stat = write( header->fp[fpndx], data, sizeof(float)*ndxcnt );
		if ( stat != sizeof(float)*ndxcnt ) {
			if ( header->debug_level & 28 ) 
			   printf( " *Err writef %i, errno=%i\n",stat,errno);
			return RADSMM_ERROR_writing;
		}
/*		if ( header->debug_level & 16 ) 
		  printf( " *** wrote out %E(%E) \n",fbuffer,*data); */
	} else
		return RADSMM_ERROR_badparam;

	return 0;
}


int RADSMM_sync( RADSMM_header_type *header )
/* this routine calls fsync */
{
	int 	i, stat;

	for ( i=0; i<header->number_of_files ; i++ ) {
		stat = fsync( header->fp[i] );
		if ( stat != 0 ) {
			if ( header->debug_level & 3)
			    printf( " fsync error %i,%i\n",stat,errno);
			return RADSMM_ERROR_writing;
		}
	}
	return 0;
}
