/* Random Access DAta Storage for MLip Mutant Version  (RADSMM) */

/* this routine opens an existing file. These routines basically fill
   in the important parts of the "header" structure. */

#include <stddef.h>
#include <unistd.h>
#include <stdio.h>
#include <strings.h>
#include <malloc.h>
#include <errno.h>
#include <string.h>
#include <sys/file.h>
#include "RADSMM.h"
#include "RADSMM_file_header.h"

#define RADSMM_VERSION 1

int RADSMM_create_file( RADSMM_header_type *header, 
			   char	*filename, int mode ) 

{
	long        i,j;
	long long   fsize, ii;
	long 	    curroff;
	double 	    dsize;
	int    	    filecnt;
	int         fp;
	ssize_t	    stat;
	double	    dbuffer[1000];
	float	    fbuffer[1000];

	RADSMM_file_header_type	fhead;

       if (header->open_flag == 'O') return RADSMM_ERROR_already_open;

/*  Step one, get the size, mostly to check */
        i = RADSMM_file_size( header, &dsize, &filecnt );
	if ( i !=0 ) {
           if ( header->debug_level & 3 )
                    printf(" *Error checking file size %li\n",i);
  	   return RADSMM_ERROR_badparam;
        }

/* Check and load the model type */
        switch ( header->model_type ) {
           case 'D' :
           case 'Q' :
	      fhead.model_type = header->model_type;
              break;
           default :
              if ( header->debug_level & 3 )
                    printf(" Bad Model_type %1c\n",header->model_type);
 	      return RADSMM_ERROR_badparam;
	}
/* Check and load the point type */
        switch ( header->marker_type ) {
           case '2' :
           case 'M' :
 	      fhead.marker_type = header->marker_type;
              break;
           default :
              if ( header->debug_level & 3 )
                    printf(" Bad marker_type %1c\n",header->marker_type);
	      return RADSMM_ERROR_badparam;
	}
/* don't bother checking the Diseq type  or LC_count*/
        fhead.use_Diseq = header->use_Diseq;
        fhead.LC_count = header->LC_count;

	if ( 	(header->element_data_type == 'F') ||
	 	(header->element_data_type == 'D') ) {
	    fhead.element_data_type = header->element_data_type;
	} else {
            if ( header->debug_level & 3 )
                    printf(" Bad element_data__type %1c\n",header->element_data_type);
	    return RADSMM_ERROR_badparam;
	}

/* Check the Header for probs, as we copy it into the structure. */
/* And calculate the offsets too.. */
	curroff = sizeof( fhead );

/* Marker (always written) */
	if ( 	(header->marker_count > 0 )&&
		(header->marker_count <= RADSMM_MAX_markers) ) {
	    fhead.marker_count = header->marker_count;
	} else {
            if ( header->debug_level & 3 )
                    printf(" Bad Marker_count %li\n",header->marker_count);
	    return RADSMM_ERROR_badparam;
	}

	fhead.marker_offset = curroff;
	header->marker_offset = fhead.pedigree_offset;
        curroff = curroff + header->marker_count*sizeof(RADSMM_marker_type);
/* Peidgree (always written) */
	if (	(header->pedigree_count > 0 ) && 
		(header->pedigree_count < RADSMM_MAX_pedigrees ) ) {
	    fhead.pedigree_count = header->pedigree_count;
	} else {
            if ( header->debug_level & 3 )
                    printf(" Bad pedigree_count %li\n",header->pedigree_count);
	    return RADSMM_ERROR_badparam;
	}
	fhead.pedigree_offset = curroff;
	header->pedigree_offset = fhead.pedigree_offset;
	curroff = curroff + header->pedigree_count*sizeof(RADSMM_pedigree_type);

/* Theta (only if marker_type=2)*/
        if ( fhead.marker_type == '2' ) {
	   if (	(header->theta_count > 0 ) && 
		(header->theta_count < RADSMM_MAX_thetas ) ) {
	      fhead.theta_count = header->theta_count;
	   } else {
              if ( header->debug_level & 3 )
                    printf(" Bad theta_count %li\n",header->theta_count);
	      return RADSMM_ERROR_badparam;
	   }
	   fhead.theta_offset = curroff;
	   header->theta_offset = fhead.theta_offset;
	   fhead.theta_matrix_type = header->theta_matrix_type;
	   curroff = curroff + header->theta_count*sizeof(RADSMM_theta_type);
        } else {
           fhead.theta_offset = 0;
           fhead.theta_matrix_type = 'D';
        }

/* Penetrance (only if model_type='D'*/
	if ( fhead.model_type == 'D' ) { 
	    if ( (header->penetrance_count > 0 ) && 
		    (header->penetrance_count < RADSMM_MAX_penetrances ) ) {
	    	fhead.penetrance_count = header->penetrance_count;
	    } else {
               if ( header->debug_level & 3 )
                    printf(" Bad penetrance_count %li\n",header->penetrance_count);
	    	return RADSMM_ERROR_badparam;
	    }
	    fhead.penetrance_offset = curroff;
	    header->penetrance_offset = fhead.penetrance_offset;
	    curroff = curroff + header->penetrance_count*sizeof(RADSMM_penetrance_type)*header->LC_count;
	} else {			 /* for models Q */
	    fhead.penetrance_count = 0;
	    fhead.penetrance_offset = 0;
	}

/* Qmodel (only if model_type='Q')*/
	if ( fhead.model_type == 'Q' ) {
	    if ( (header->qmodel_count > 0 ) && 
		    (header->qmodel_count < RADSMM_MAX_qmodels ) ) {
	    	fhead.qmodel_count = header->qmodel_count;
	    } else {
              if ( header->debug_level & 3 )
                    printf(" Bad qmodel_count %li\n",header->qmodel_count);
	    	return RADSMM_ERROR_badparam;
	    }
	    fhead.qmodel_offset = curroff;
	    header->qmodel_offset = fhead.qmodel_offset;
	    curroff = curroff + header->qmodel_count*sizeof(RADSMM_quantitative_type);
	} else {
	    fhead.qmodel_count = 0;
	    fhead.qmodel_offset = 0;
	}

/* dieseq (only if Use_diseq!='N'*/
	if ( fhead.use_Diseq !=  'N' ) {
	    if ( (header->diseq_count > 0 ) && 
		    (header->diseq_count < RADSMM_MAX_diseq_params ) ) {
	    	fhead.diseq_count = header->diseq_count;
	    } else {
                if ( header->debug_level & 3 )
                    printf(" Bad diseq_count %li\n",header->diseq_count);
	    	return RADSMM_ERROR_badparam;
	    }
           fhead.diseq_offset = curroff;
           header->diseq_offset = fhead.diseq_offset;
           curroff = curroff + header->diseq_count*sizeof(RADSMM_diseq_type);
	} else {
	    fhead.diseq_count = 0;
	    fhead.diseq_offset = 0;
	}

/* geneFreq (only if model_type='D') */
	if ( fhead.model_type ==  'D' ) {
 	   if (	(header->geneFreq_count > 0 ) && 
		(header->geneFreq_count < RADSMM_MAX_geneFreqs ) ) {
	      fhead.geneFreq_count = header->geneFreq_count;
	   } else {
	      return RADSMM_ERROR_badparam;
	   }
	   fhead.geneFreq_offset = curroff;
	   header->geneFreq_offset = fhead.geneFreq_offset;
	   curroff = curroff + header->geneFreq_count * sizeof(RADSMM_geneFreq_type);
        } else {
           fhead.geneFreq_count = 0;
           fhead.geneFreq_offset = 0;
        }

/* markerlabel */
	fhead.markerlabel_size = header->markerlabel_size;
	fhead.markerlabel_offset = curroff;
	curroff = curroff + (header->marker_count * header->markerlabel_size);
	/*header->markerlabel_offset = fhead.markerlabel_offset; */

/* pedigreelabel */
	fhead.pedigreelabel_size = header->pedigreelabel_size;
	fhead.pedigreelabel_offset = curroff;
	curroff = curroff + (header->pedigree_count * header->pedigreelabel_size);
	/*header->pedigreelabel_offset = fhead.pedigreelabel_offset;*/
/* */
	if ( header->debug_level & 16 )
	    printf(" *** Markerlabel cnt:%li, size:%li Pedlabel cnt:%li, size:%li\n",
			header->marker_count, header->markerlabel_size, 
			header->pedigree_count, header->pedigreelabel_size );
	fhead.start_of_data = curroff;

	header->start_of_data = fhead.start_of_data;

	fhead.number_of_files = header->number_of_files;
	fhead.chunks_per_file =	header->chunks_per_file;

	if ( header->debug_level & 16 )
	    printf(" *** Offset values ped:%li, theta:%li, penet:%li, GF:%li, Markerlabel:%li, Pedlabel:%li, SOD:%li\n",
				fhead.pedigree_offset, fhead.theta_offset,
				fhead.penetrance_offset, fhead.geneFreq_offset,
				fhead.markerlabel_offset, fhead.pedigreelabel_offset,
				fhead.start_of_data );
 
	/* version number */
	header->version = RADSMM_VERSION;
	sprintf( fhead.cookie, RADSMM_cookie );
	fhead.version = header->version;
	fhead.subversion = 0;
	fhead.ordering = header->ordering;

	memcpy( fhead.date_string, header->date_string, 17);
	memcpy( fhead.description, header->description, 64);

/* Check for file existance, but opening the file as excl */

	/* Ok, Create the lock the file */
	fp = open( filename, O_RDWR|O_CREAT|O_EXCL/*|O_LARGEFILE*/, mode );
	if ( fp == -1 ) {
		if ( header->debug_level & 1)
			printf(" open file error %i\n",errno);
		return RADSMM_ERROR_fileopen;
	}

	if ( lockf( fp, F_LOCK, 0) == -1 ) {
		if ( (header->debug_level & 1) != 0 )
			printf(" locking error %i",errno);
		close ( fp );
		return RADSMM_ERROR_locking;
	}

	if ( header->debug_level & 16 ) printf(" *** File Created and locked.\n");

	/* store fp in the header */
	header->fp[0] = fp;

/* Write the header */
	stat = write( fp, &fhead, sizeof(fhead) );
	if ( stat != sizeof(fhead) )  {
/*		if ( header->debug_level & 1)  */
			printf("* write header error %i:%i:%s",stat,errno,strerror(errno));
		lockf( fp, F_ULOCK, 0 );
		close(fp);
		return RADSMM_ERROR_writing;
	}

/* Now write the indexes and check the file offsets too */
	/* marker */
	if ( lseek(fp, 0L, SEEK_CUR) != fhead.marker_offset) {
		printf( " Error: marker offset %li %li",fhead.marker_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	}
	stat =write( fp, header->marker_list, 
                           header->marker_count*sizeof(RADSMM_marker_type));
	if ( stat <=0 ) {
		if ( header->debug_level & 1)
			printf(" Marker index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf ( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
	}

	/* Pedigree */
	if ( lseek(fp, 0L, SEEK_CUR) != fhead.pedigree_offset) {
		printf( " Error: pedigree offset %li %li",fhead.pedigree_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0 );
		close(fp);
		return RADSMM_ERROR_lseek;
	}
	stat = write( fp, header->pedigree_list, 
                           header->pedigree_count*sizeof(RADSMM_pedigree_type) );
	if ( stat <=0 ) {
/*		if ( header->debug_level & 1) */
			printf("* pedigree index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf ( fp, F_ULOCK, 0 );
		close(fp);
		return RADSMM_ERROR_writing;
	}

	/* theta */
        if ( header->marker_type == '2' ) {
 	   if ( lseek(fp, 0L, SEEK_CUR) != fhead.theta_offset) {
		printf(" Error: theta offset %li %li",fhead.theta_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	   }
	   stat = write( fp, header->theta_list, 
                           header->theta_count*sizeof(RADSMM_theta_type));
	   if ( stat <=0 ) {
		if ( header->debug_level & 1)
			printf(" theta index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
           }
	}

	/* penetrance */
	if ( fhead.penetrance_offset > 0 ) {
	    if ( lseek(fp, 0L, SEEK_CUR) != fhead.penetrance_offset) {
		printf(" Error: penetrance offset %li %li",fhead.penetrance_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	    }
            for ( i=0 ; i<header->LC_count ; i++ ) {
	       stat = write( fp, header->penetrance_list[i], header->penetrance_count*sizeof(RADSMM_penetrance_type) ) ;
            }
	    if ( stat <=0 ) {
		if ( header->debug_level & 1) 
			printf(" penet index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
	    }
	}

	/* qmodel */
	if ( fhead.qmodel_offset > 0 ) {
	    if ( lseek(fp, 0L, SEEK_CUR) != fhead.qmodel_offset) {
		printf(" Error: qmodel offset %li %li",fhead.qmodel_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	    }
	    stat = write( fp, header->qmodel_list, header->qmodel_count*sizeof(RADSMM_quantitative_type)) ;
	    if ( stat <=0 ) {
		if ( header->debug_level & 1) 
			printf(" mean index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
	    }
	}


	/* diseq */
        if ( header->use_Diseq != 'N' ) {
 	   if ( lseek(fp, 0L, SEEK_CUR) != fhead.diseq_offset) {
		printf(" Error: Diseq offset %li %li",fhead.penetrance_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	   }
	   stat = write( fp, header->diseq_list, header->diseq_count*sizeof(RADSMM_diseq_type)) ;
	   if ( stat <=0 ) {
		if ( header->debug_level & 1) 
			printf(" diseq index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
           }
	}


	/* geneFreq */
	if ( fhead.geneFreq_offset > 0 ) {
	   if ( lseek(fp, 0L, SEEK_CUR) != fhead.geneFreq_offset) {
		printf( " Error: geneFreq offset %li %li",fhead.geneFreq_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	   }
	   stat =write( fp, header->geneFreq_list, 
                           header->geneFreq_count*sizeof(RADSMM_geneFreq_type));
	   if ( stat <=0 ) {
		if ( header->debug_level & 1)
			printf(" genefreq index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf ( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
	   }
        }

	/* marker labels */
	if ( header->markerlabel_size > 0 ) {
	    if ( lseek(fp, 0L, SEEK_CUR) != fhead.markerlabel_offset) {
		printf( " Error: markerlabel offset %li %li",fhead.markerlabel_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	    }
	    stat =write( fp, header->markerlabel_list, 
                           header->marker_count*header->markerlabel_size);
	    if ( stat <=0 ) {
		if ( header->debug_level & 1)
			printf(" markerlabel index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf ( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
	    }
	}

	/* pedigree labels */
	if ( header->pedigreelabel_size > 0 ) {
	    if ( lseek(fp, 0L, SEEK_CUR) != fhead.pedigreelabel_offset) {
		printf( " Error: pedigreelabel offset %li %li",fhead.pedigreelabel_offset, lseek(fp, 0, SEEK_CUR));
		lockf( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_lseek;
	    }
	    stat =write( fp, header->pedigreelabel_list, 
                           header->pedigree_count*header->pedigreelabel_size);
	    if ( stat <=0 ) {
		if ( header->debug_level & 3)
			printf(" pedigreelabel index write error %i:%i:%s",stat,errno,strerror(errno));
		lockf ( fp, F_ULOCK, 0L );
		close(fp);
		return RADSMM_ERROR_writing;
	    }
	}

	if ( header->debug_level & 16 ) 
		printf("*** File header written\n" );


	/* Clear out the buffers we are going to write to the file(s) */
	for ( i=0; i<1000; i++) {
		dbuffer[i] = (double) RADSMM_EMPTY;
		fbuffer[i] = (float) RADSMM_EMPTY;
	}

	if ( header->number_of_files == 0 ) {
	    /* lets check the final offset */
	    if ( lseek(fp, 0L, SEEK_CUR) != header->start_of_data ) {
		if ( header->debug_level & 3)
		   printf( " Error: startofdata offset %li %li",fhead.penetrance_offset, lseek(fp, 0, SEEK_CUR));
		return RADSMM_ERROR_lseek;;
	    }
	
	    /* Now to fill the file with blanks */
	    /* Note that it writes markers in group. the old way, a faster way */
	    fsize = header->marker_count * header->pedigree_count * header->theta_count 
		* header->penetrance_count * header->geneFreq_count * header->qmodel_count * header->diseq_count;
	    if ( header->theta_matrix_type == 'G') 
		fsize = fsize * header->theta_count;

	    for ( ii=0; ii<fsize; ii+=1000 ) {
	      if ( ii+1000 > fsize ) {
		  j=fsize-ii+1;
	      } else {
		j = 1000;
	      }
	      if ( header->element_data_type == 'D' ) {
		stat = write( fp, dbuffer, j*sizeof(double) );
		if ( stat < 0 ) {
		    if ( header->debug_level & 3)
			printf( " Error: writing empty %i %lli %li\n",stat, ii,j);
		    lockf( fp, F_ULOCK, 0 );
		    close(fp);
		    return RADSMM_ERROR_writing;
		}
	      } else if ( header->element_data_type == 'F' ) {
		stat = write( fp, fbuffer, j*sizeof(float) );
		if ( stat < 0 ) {
		    if ( header->debug_level & 3)
		       printf( " Error: writing empty %i %lli %li\n",stat, ii, j);
		    lockf( fp, F_ULOCK, 0 );
		    close(fp);
		    return RADSMM_ERROR_writing;
		}
	      } else {
		lockf( fp, F_ULOCK, 0 );
		close(fp);
		return RADSMM_ERROR_internal;
	      }
	    }

	} else {  /* more than 1 file */
	    dbuffer[0] = (double) RADSMM_EMPTY;
	    fbuffer[0] = (float) RADSMM_EMPTY;
	    for ( i=1 ; i<=header->number_of_files ; i++ ) {
		if ( header->fp[i] < 0 ) {
		    if ( header->debug_level & 3 )printf( " Error: opening extra file \"%s\",%i\n",filename,stat);
		    lockf( fp, F_ULOCK, 0 );
		    close(fp);
		    return RADSMM_ERROR_fileopen;
		}

		for ( ii=0 ; ii<header->chunks_per_file ; ii+=1000 ) {
 		    if ( ii+1000 > header->chunks_per_file ) {
			j = header->chunks_per_file - ii + 1;
		    } else {
			j = 1000;
		    }
		    if ( header->element_data_type == 'D') {
			stat = write( header->fp[i], dbuffer, j*sizeof(double) );
			if ( stat < 0 ) {
			    if ( header->debug_level & 3 )
				printf( " Error: writing empty %i f=%li, cnt=%lli\n",stat, i, ii);
			    lockf( fp, F_ULOCK, 0 );
			    close(fp);
			    return RADSMM_ERROR_writing;
			}
		    } else {
			stat = write( header->fp[i], fbuffer, j*sizeof(float) );
			if ( stat < 0 ) {
			    if ( header->debug_level & 3 )
				printf( " Error: writing empty %i f=%li, cnt=%lli\n",stat, i, ii);
			    lockf( fp, F_ULOCK, 0 );
			    close(fp);
			    return RADSMM_ERROR_writing;
			}
		    }
		}
	    }
	}
	if ( header->debug_level & 16 ) 
		printf("*** successful file creation\n" );
 
	header->open_flag = 'O';
	return RADSMM_ERROR_success;

}
