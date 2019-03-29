/* Random Access DAta Storage for MLip Mutant Version (RADSMM) */


/* this routine opens an existing file. These routines basically fill
   in the important parts of the "header" structure. */


#include <stddef.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <malloc.h>
#include "RADSMM.h"
#include "RADSMM_file_header.h"
#include <sys/file.h>
#include <errno.h>

#define RADSMM_VERSION 1

int RADSMM_close_file( RADSMM_header_type *header )
  {
	int	i,rc;
	int     stat = 0;

	if ( header->open_flag != 'O' ) return RADSMM_ERROR_not_open;

	for ( i=0 ; i<=header->number_of_files ; i++ ) {
		lockf( header->fp[i], F_ULOCK, 0L );   /* note that we ignore error codes from the unlock */
		stat = close(header->fp[i]);
		if ( stat != 0 ) { 
		   rc = errno;
		   if ( header->debug_level & 3 ) 
		      printf(" Error closing file %i, err=%i",header->fp[i],errno);
		}
	}
	header->open_flag = 'N';
	return stat;
}  

int RADSMM_open_file( RADSMM_header_type *header, 
			   char	*filename, char open_mode ) 
{
	int    i;
	int    fp;
	ssize_t	stat;
/*	int	version;  */
	RADSMM_file_header_type	fhead;
	int 	mode;


	if ( header->open_flag == 'O' ) return RADSMM_ERROR_already_open;

/* opening the file  */
	if ( (open_mode == 'W') || (open_mode == 'w') ) 
		mode = O_RDWR;
	else
		mode = O_RDONLY;

	/* Ok, open and lock the file */
	fp = open( filename, mode, 0 );
	if ( fp == -1 ) {
		if ( header->debug_level & 3 )
		    printf(" Error opening file %i\n",errno );
		return RADSMM_ERROR_fileopen;
	}

	if ( mode != O_RDONLY ) {
	    if ( lockf( fp, F_LOCK ,0 ) == -1 ) {
		if ( header->debug_level & 3 )
		    printf(" Error locking file %i, %i\n",errno,fp );
		close ( fp );
		return RADSMM_ERROR_locking;
	    }
	}
	/* store fp in the header */
	header->fp[0] = fp;
	header->number_of_files = 1;

/* Read the header */
	stat = read( fp, &fhead, sizeof(fhead) );
	if ( stat<=0 )  {
		if ( header->debug_level & 3 )
		    printf(" Error in Open, reading header %i,%i\n",stat,errno); 
		lockf( fp, F_ULOCK, 0 );
		close(fp);
		return RADSMM_ERROR_reading;
	}


/* Check the Cookie */
	if ( strncmp( fhead.cookie, RADSMM_cookie, 4) == 0 ) {
		/* everything is ok, do nothing */
	} else {
		if ( header->debug_level & 3 ) 
			printf( " File is not a RADSMM file. \"%.4s\"\n",fhead.cookie);
		return RADSMM_ERROR_file_header;
	}


/* Check the version */
	header->version = fhead.version;
	header->subversion = fhead.subversion;
	switch( header->version ) {
	    case 1:
		/* supported, do nothing */
		break;
	    default:
		if ( header->debug_level & 3 ) 
		    printf(" File Version not supported.:%i",header->version);
		return RADSMM_ERROR_file_header;
	}


/* open other files if need */
	if ((fhead.number_of_files<0) || (fhead.number_of_files>MAX_DATA_FILES)) {
	    if ( header->debug_level & 3 ) 
		printf( " ERROR: header.Number of files = %i.\n", fhead.number_of_files );
	    return RADSMM_ERROR_file_header;
	}
	header->number_of_files = fhead.number_of_files;

	if ( fhead.chunks_per_file <= 0 ) {
	    if ( header->debug_level & 3 ) 
		printf( " ERROR: header.chunks per file = %li.\n", fhead.chunks_per_file );
	    return RADSMM_ERROR_file_header;
	}
	header->chunks_per_file = fhead.chunks_per_file;

	if ( header->number_of_files >= 1 ) {
		if ( header->debug_level & 16 )
		   printf(" ***Opening %i extra data files.\n",header->number_of_files); 
		for ( i=1 ; i <= header->number_of_files ; i ++ ) {
			header->fp[i] = RADSMM_open_others( filename, i, mode, 0 );
			if ( header->fp[i] < 0 ) {
			    if ( header->debug_level & 3 )
				   printf("Error opening data file %i: %i,%i\n",i,stat,errno);
			    return RADSMM_ERROR_fileopen;
			}
			/* once the files is open, we _DONT_ check it */
		}
	} 


/* Check the fileHeader for probs, as we copy it into the structure. */

	header->chunks_per_file = fhead.chunks_per_file;
	if ( fhead.element_data_type == 'F') {
	    	header->chunk_size = sizeof( float );
		header->element_data_type = 'F';
	} else if (fhead.element_data_type == 'D') {
	    	header->chunk_size = sizeof( double );
		header->element_data_type = 'D';
	} else {
	    if ( header->debug_level & 3) 
		printf(" File Header Problem: Bad data type %.1s\n",
			&fhead.element_data_type);
	    return RADSMM_ERROR_file_header;
	}
        switch ( fhead.model_type ) {
           case 'Q':
           case 'D':
              header->model_type = fhead.model_type;
              break;
           default:
	      if ( header->debug_level & 3) 
	         printf(" File Header Problem: Bad model type %.1s\n",
			&fhead.model_type);
	      return RADSMM_ERROR_file_header;
        }
        switch ( fhead.marker_type ) {
           case '2':
           case 'M':
              header->marker_type = fhead.marker_type;
              break;
           default:
	      if ( header->debug_level & 3) 
	         printf(" File Header Problem: Bad marker type %.1s\n",
			&fhead.element_data_type);
	      return RADSMM_ERROR_file_header;
        }
        /* no real checks on use_Diseq */
        header->use_Diseq = fhead.use_Diseq;

   /* lets do marker stuff */
	if ( 	(fhead.marker_count > 0 )&&
		(fhead.marker_count <= RADSMM_MAX_markers) ) {
	    header->marker_count = fhead.marker_count;
	} else {
	    if ( header->debug_level & 3) 
		printf(" File Header Problem: Bad elements per Marker %li\n",
			fhead.marker_count);
	    return RADSMM_ERROR_file_header;
	}
/*	header->marker_offset = fhead.marker_offset; */
   
        if ((header->marker_list = malloc(header->marker_count*
					sizeof(RADSMM_marker_type))) == NULL ) {
	    if ( header->debug_level & 3 )
		printf(" Error allocating memory for Markers %i",errno);
            return RADSMM_ERROR_malloc;
	}

	stat = lseek(header->fp[0], fhead.marker_offset, SEEK_SET);
        stat = read( header->fp[0], header->marker_list, 
				header->marker_count*sizeof(RADSMM_marker_type) );
	if ( stat <= 0 ) {
	    if ( header->debug_level & 3 )
		printf(" Error reading file for marker list %i:%i",stat,errno);
            return RADSMM_ERROR_reading;
        }


   /* lets do pedigree stuff */
	if (	(fhead.pedigree_count > 0 ) && 
		(fhead.pedigree_count < RADSMM_MAX_pedigrees ) ) {
	    header->pedigree_count = fhead.pedigree_count;
	} else {
	    if ( header->debug_level & 3) 
		printf(" File Header Problem: Bad pedigree cnt %li\n",
				fhead.pedigree_count);
	    return RADSMM_ERROR_file_header;
	}
/*	header->pedigree_offset = fhead.pedigree_offset; */
   
        if ((header->pedigree_list = malloc(header->pedigree_count*
					sizeof(RADSMM_pedigree_type))) == NULL ) {
	    if ( header->debug_level & 3 )
		printf(" Error allocating memory for pedigree %i",errno);
            return RADSMM_ERROR_malloc;
	}

	stat = lseek(header->fp[0], fhead.pedigree_offset, SEEK_SET);
        stat = read( header->fp[0], header->pedigree_list, 
				header->pedigree_count*sizeof(RADSMM_pedigree_type) );
	if ( stat <= 0 ) {
	    if ( header->debug_level & 3 )
		printf(" Error reading file for pedigreelist %i:%i",stat,errno);
            return RADSMM_ERROR_reading;
	}
		
   /* lets do theta stuff */
        if ( header->marker_type != '2' ) {
           header->theta_count = 1;
           header->theta_matrix_type = 'D';
        } else {
	   if ( (fhead.theta_count > 0 ) && 
		(fhead.theta_count < RADSMM_MAX_thetas ) ) {
	      header->theta_count = fhead.theta_count;
	   } else {
	      printf(" File Header Problem: Bad theta cnt %li\n",
				fhead.theta_count);
	      return RADSMM_ERROR_file_header;
	   }
/*	   header->theta_offset = fhead.theta_offset; */
   
           if ((header->theta_list = malloc(header->theta_count*
					sizeof(RADSMM_theta_type))) == NULL ) {
	      if ( header->debug_level & 3 )
		         printf(" Error allocating memory for theta %i",errno);
              return RADSMM_ERROR_malloc;
	   }

	   stat = lseek(header->fp[0], fhead.theta_offset, SEEK_SET);
           stat = read( header->fp[0], header->theta_list, 
				header->theta_count*sizeof(RADSMM_theta_type) );
	   if ( stat <= 0 ) {
	      if ( header->debug_level & 3 )
		printf(" Error reading file for theta list %i:%i",stat,errno);
              return RADSMM_ERROR_reading;
	   }

	   if ( (fhead.theta_matrix_type == 'D') || (fhead.theta_matrix_type == 'G') ) {
              header->theta_matrix_type = fhead.theta_matrix_type;
	   } else {
	      if ( header->debug_level & 3 ) 
		    printf(" Error in theta Matrix type \"%s\",%i\n",&fhead.theta_matrix_type,fhead.theta_matrix_type);
	      return RADSMM_ERROR_file_header;
	   }
        }

   /* lets do penetrence stuff   (Now with LC!!) */
        if ( header-> model_type != 'D' ) {
           header->penetrance_count = 1;
           header->LC_count = 0;
        } else {
	   if (	(fhead.penetrance_count > 0 ) && 
		(fhead.penetrance_count < RADSMM_MAX_penetrances ) ) {
	      header->penetrance_count = fhead.penetrance_count;
	   } else {
              if ( header->debug_level & 3) 
		 printf(" File Header Problem: Bad penetrance cnt %li\n",
				fhead.penetrance_count);
	      return RADSMM_ERROR_file_header;
	   }
/*	   header->penetrance_offset = fhead.penetrance_offset; */

	   if ( (fhead.LC_count > 0 ) && 
		    (fhead.LC_count < RADSMM_MAX_liability_classes ) ) {
	      header->LC_count = fhead.LC_count;
	   } else {
              if ( header->debug_level & 3) 
		 printf(" File Header Problem: Bad LC cnt %li\n",
				fhead.LC_count);
	         return RADSMM_ERROR_file_header;
	   }
   
	   stat = lseek(header->fp[0], fhead.penetrance_offset, SEEK_SET);


           for ( i=0 ; i< header->LC_count ; i++ ) {
              if ((header->penetrance_list[i] = malloc(header->penetrance_count*
			header->LC_count*sizeof(RADSMM_penetrance_type))) == NULL ) {
	         if ( header->debug_level & 3 )
		    printf(" Error allocating memory for penetrance %i",errno);
                 return RADSMM_ERROR_malloc;
	      }
              stat = read( header->fp[0], header->penetrance_list[i], 
			header->penetrance_count*sizeof(RADSMM_penetrance_type) );
	      if ( stat <= 0 ) {
	         if ( header->debug_level & 3 )
	            printf(" Error reading file for penetrancelist for lc= %i  %i:%i\n",i, stat,errno);
                 return RADSMM_ERROR_reading;
              }
           }
        }

   /* lets do qmodel stuff */
	if ( header->model_type != 'Q' ) {
           header->qmodel_count = 1;
        } else {
           if ( (fhead.qmodel_count > 0 ) && 
		      (fhead.qmodel_count < RADSMM_MAX_qmodels ) ) {
	      header->qmodel_count = fhead.qmodel_count;
	   } else {
              if ( header->debug_level & 3) 
	           printf(" File Header Problem: Bad qmodel cnt %li\n",
				fhead.qmodel_count);
	      return RADSMM_ERROR_file_header;
	   }
	   header->qmodel_offset = fhead.qmodel_offset;
   
           if ((header->qmodel_list = malloc(header->qmodel_count*
					sizeof(RADSMM_quantitative_type))) == NULL ) {
	      if ( header->debug_level & 3 )
	             printf(" Error allocating memory for qmodel %i",errno);
              return RADSMM_ERROR_malloc;
	   }

	   stat = lseek(header->fp[0], fhead.qmodel_offset, SEEK_SET);
           stat = read( header->fp[0], header->qmodel_list, 
				header->qmodel_count*sizeof(RADSMM_quantitative_type) );
	   if ( stat <= 0 ) {
	      if ( header->debug_level & 3 )
		printf(" Error reading file for qmodel list %i:%i",stat,errno);
              return RADSMM_ERROR_reading;
           }
        }

   /* lets do diseq stuff */
        if ( header->use_Diseq == 'N' ) {
           header->diseq_count = 1;
        } else {
	   if ( (fhead.diseq_count > 0 ) && 
		(fhead.diseq_count < RADSMM_MAX_diseq_params ) ) {
	      header->diseq_count = fhead.diseq_count;
	   } else {
              if ( header->debug_level & 3) 
		printf(" File Header Problem: Bad diseq cnt %li\n",
				fhead.diseq_count);
	      return RADSMM_ERROR_file_header;
	   }
/*	   header->diseq_offset = fhead.diseq_offset; */
   
           if ((header->diseq_list = malloc(header->diseq_count*
					sizeof(RADSMM_diseq_type))) == NULL ) {
	      if ( header->debug_level & 3 )
		printf(" Error allocating memory for diseq %i",errno);
              return RADSMM_ERROR_malloc;
	   }

	   stat = lseek(header->fp[0], fhead.diseq_offset, SEEK_SET);
           stat = read( header->fp[0], header->diseq_list, 
				header->diseq_count*sizeof(RADSMM_diseq_type) );
	   if ( stat <= 0 ) {
	      if ( header->debug_level & 3 )
		printf(" Error reading file for diseq list %i:%i",stat,errno);
              return RADSMM_ERROR_reading;
           }
        }

   /* lets do geneFreq stuff */
        if ( header->model_type != 'D' ) {
           header->geneFreq_count = 1;
        } else {
	   if ( (fhead.geneFreq_count > 0 ) && 
		(fhead.geneFreq_count < RADSMM_MAX_geneFreqs ) ) {
	      header->geneFreq_count = fhead.geneFreq_count;
	   } else {
              if ( header->debug_level & 3) 
		printf(" File Header Problem: Bad geneFreq cnt %li\n",
				fhead.geneFreq_count);
	      return RADSMM_ERROR_file_header;
	   }
/*	   header->geneFreq_offset = fhead.geneFreq_offset; */
   
           if ((header->geneFreq_list = malloc(header->geneFreq_count*
					sizeof(RADSMM_geneFreq_type))) == NULL ) {
	      if ( header->debug_level & 3 )
		printf(" Error allocating memory for geneFreq %i",errno);
              return RADSMM_ERROR_malloc;
	   }

	   stat = lseek(header->fp[0], fhead.geneFreq_offset, SEEK_SET);
           stat = read( header->fp[0], header->geneFreq_list, 
				header->geneFreq_count*sizeof(RADSMM_geneFreq_type) );
	   if ( stat <= 0 ) {
	      if ( header->debug_level & 3 )
		printf(" Error reading file for geneFreqlist %i:%i",stat,errno);
              return RADSMM_ERROR_reading;
           }
        }


   /* lets do the Marker labels */
	if ( fhead.markerlabel_size < 0 )  {
	    if ( header->debug_level & 3 ) 
		printf(" File Header Problem: Bad markerlabel size %li\n",
				fhead.markerlabel_size );
	    return RADSMM_ERROR_file_header;
	} else if ( fhead.markerlabel_size == 0 ) {
		/* no data, do nothing */
		header->markerlabel_list = NULL;
	} else {
	    header->markerlabel_size = fhead.markerlabel_size;
	    if ((header->markerlabel_list = calloc(header->marker_count*
			header->markerlabel_size, sizeof(char))) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" Error allocating memory for markerlabel %i, (%li,%li)\n",
			errno, header->marker_count, header->markerlabel_size);
        	return RADSMM_ERROR_malloc;
	    }

	    stat = lseek(header->fp[0], fhead.markerlabel_offset, SEEK_SET);
	    stat = read( header->fp[0], header->markerlabel_list, 
			header->marker_count*header->markerlabel_size );
	    if ( stat <= 0 ) {
		if ( header->debug_level & 3 )
		    printf(" Error reading file for markerlabel %i:%i",stat,errno);
        	return RADSMM_ERROR_reading;
	    }
	}

   /* Lets do the Pedigree labels */
	if ( fhead.pedigreelabel_size < 0 )  {
	    if ( header->debug_level & 3 ) 
		printf(" File Header Problem: Bad pedigreelabel size %li\n",
				fhead.pedigreelabel_size);
	    return RADSMM_ERROR_file_header;
	} else if ( (fhead.pedigreelabel_size == 0) ) {
		/* no data, do nothing */
		header->pedigreelabel_list = NULL;
	} else {
	    header->pedigreelabel_size = fhead.pedigreelabel_size;
	    if ((header->pedigreelabel_list = calloc(header->pedigree_count*
			header->pedigreelabel_size, sizeof(char))) == NULL) {
		if ( header->debug_level & 3 )
		    printf(" Error allocating memory for pedigreelabel %i",errno);
        	return RADSMM_ERROR_malloc;
	    }

	    stat = lseek(header->fp[0], fhead.pedigreelabel_offset, SEEK_SET);
	    stat = read( header->fp[0], header->pedigreelabel_list, 
			header->pedigree_count*header->pedigreelabel_size );
	    if ( stat <= 0 ) {
		if ( header->debug_level & 3 )
		    printf(" Error reading file for pedigreelabel %i:%i\n",stat,errno);
        	return RADSMM_ERROR_reading;
	    }
	}

	/* grab some more stuff */
	header->start_of_data = fhead.start_of_data;
/*	if ( header->debug_level &16 ) 
		printf( "start of data = %i",header->start_of_data); */

	/* set some other stuff */
	header->open_flag = 'O';
	if ( fhead.version == 1 ) {
		header->ordering = 'A';
	} else {
		header->ordering = fhead.ordering;
	}

	memcpy( header->date_string, fhead.date_string, 17);
	memcpy( header->description, fhead.description, RADSMM_descr_length);

	header->open_flag='O';
	return 0;

}


int RADSMM_open_others( char filename[], int number, int flag, int open_mode  )
/* THis routine opens the other files in a multi file dataset */

{
	int	fp;

	char suffix[3];
	char  newname[1024];

	sprintf( suffix, "_%2.2i", number);

	strcpy( newname, filename );
	strncat( newname, suffix, 3 );

/*	if ( header->debug_level & 16 ) */
	    printf( " *** new file name \"%s\"\n", newname );
 
	/* Ok, open and lock the file */
	fp = open( newname, flag, open_mode );
	if ( fp == -1 ) {
/*		if ( header->debug_level & 3 ) */
		    printf(" Error opening file %i\n",errno );
		return RADSMM_ERROR_fileopen;
	}

	if ( flag != O_RDONLY ) {
	    if ( lockf( fp, F_LOCK ,0 ) == -1 ) {
/*		if ( header->debug_level & 3 ) */
		    printf(" Error locking file %i, %i\n",errno,fp );
		close ( fp );
		return RADSMM_ERROR_locking;
	    }
	}
	return fp;
}
