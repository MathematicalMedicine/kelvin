 /* Random Access DAta Storage for MLip (Multiple Model)  (RADSMM) */


/* this routine sets up the header so one can create a new file.

   */


#include <stddef.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <malloc.h>
#include "RADSMM.h"


int RADSMM_setup_init( RADSMM_header_type *header, int debug_level )
{
   int  i;

	header->marker_count = 0;
	header->element_data_type = ' ';
	header->pedigree_count = 0;
	header->theta_count = 0;
	header->theta_matrix_type = ' ';
	header->LC_count = 1;
	header->penetrance_count = 0;
	header->qmodel_count = 0;
	header->diseq_count = 0;
	header->geneFreq_count = 0;
	
	header->number_of_files = 0;
	header->version = 1;
	header->subversion = 1;

	header->open_flag = 'N';
	header->start_of_data = 0;

	header->debug_level = debug_level;

	header->description[0] = '\0';
	header->markerlabel_size = 0;
	header->pedigreelabel_size = 0;

	header->checksum_type = 0;

	header->ordering = 'A';
	header->LC_count = 1;
        header->marker_type = '2';
	header->model_type = 'D';
        header->use_Diseq = 'N';

/*	for ( i = 0 ; i < RADSMM_MAX_diseq_params ; i++ ) {
           header->diseq_list[i] = NULL;
        } */
        header->diseq_list = NULL;

	for ( i = 0 ; i < RADSMM_MAX_liability_classes ; i++ ) {
           header->penetrance_list[i] = NULL;
        } 

	if ( header->debug_level & 16 ) {
		printf( " ***Debug level set to %i\n", header->debug_level );
	}
	return 0;
}

int RADSMM_setup_type( RADSMM_header_type *header, char Point_Type, 
                                              char model_type, char use_diseq )
{

   switch (Point_Type) {
      case '2':
         header->marker_type = '2';
         break;
      case 'M':
      case 'm':
         header->marker_type = 'M';
         header->theta_count = 1;
         break;
      default:
         if ( header->debug_level & 3 )
            printf(" *Point_type  for setup_type bad %c\n",Point_Type);
         return RADSMM_ERROR_badparam;
   }


   switch (model_type) {
     case 'D':
     case 'd':
        header->model_type = 'D';
        header->qmodel_count = 1;
        break;
     case 'Q':
     case 'q':
        header->model_type = 'Q';
        header->geneFreq_count = 1;
        header->penetrance_count = 1;
        break;
     default:
        if ( header->debug_level & 3 )
           printf(" *model_type for setup_type bad %c\n",Point_Type);
        return RADSMM_ERROR_badparam;
   }

   switch (use_diseq) {
     case 'N':
     case 'n':
        header->use_Diseq = 'N';
        header->diseq_count = 1;
        break;
     default:
        header->use_Diseq = use_diseq;
   }

   return 0;
}

                 
int RADSMM_setup_marker( RADSMM_header_type *header, 
			   RADSMM_marker_type *list, 
			   long count )
{
	int   i;
	RADSMM_marker_type   *p;

	if ( (count <= 0) || (count > RADSMM_MAX_markers) ) {
		if ( header->debug_level & 3 )
		    printf(" *Marker Count for RADSMM_setup_marker bad %li",count);
		return RADSMM_ERROR_badparam;
	}

	if ((header->marker_list = malloc(count*sizeof(RADSMM_marker_type))) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" *Malloc error for marker list\n");
		return RADSMM_ERROR_malloc;
	}

	header->marker_count = count;

	if ( list == NULL ) {
		if ( header->debug_level &16 ) 
		    fprintf(stderr," ***No list passed to setup_marker, assuming sequential\n");
		for ( i=0,p=header->marker_list; i<count ; i++,p++ )
			*p = i;

	} else {
 	    memcpy( header->marker_list, list, count*sizeof(RADSMM_marker_type) );
	}

	return 0;
}
int RADSMM_setup_pedigree( RADSMM_header_type *header, 
			   int *list, 
			   long count )
{
	int   i;
	int   *p;

	if ( (count <= 0) || (count > RADSMM_MAX_pedigrees) ) {
		if ( header->debug_level & 3 )
		    printf(" *Count for setup_pedigree bad %li\n",count);
		return RADSMM_ERROR_badparam;
	}

	if ((header->pedigree_list = malloc(count*sizeof(RADSMM_pedigree_type))) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" *Malloc error for pedigree list\n");
		return RADSMM_ERROR_malloc;
	}

	header->pedigree_count = count;

	if ( list == NULL ) {
		if ( header->debug_level &16 ) 
		    fprintf(stderr," ***No list passed to setup_ped, assuming sequential\n");
		for ( i=0,p=header->pedigree_list; i<count ; i++,p++ )
			*p = i;

	} else {
 	    memcpy( header->pedigree_list, list, count*sizeof(RADSMM_pedigree_type) );
	}

	return 0;
}

int RADSMM_setup_theta( RADSMM_header_type *header, 
			   double list[], 
			   long   count,
			   char thetatype )
{
	if ( (count <= 0) || (count > RADSMM_MAX_thetas) ) {
		if ( header->debug_level & 3 )
		    printf(" *Count for setup_theta bad %li",count);
		return RADSMM_ERROR_badparam;
	}


	if ( (thetatype == 'D') || (thetatype == 'd') )
		header->theta_matrix_type = 'D';
	else if ( (thetatype == 'G') || (thetatype == 'g') )
		header->theta_matrix_type = 'G';
	else {
		if ( header->debug_level & 3 )
		    printf(" *matrix type for setup_theta bad \"%1s\"\n",
					&thetatype);
		return RADSMM_ERROR_badparam;
	}

	if ((header->theta_list = malloc(count*sizeof(RADSMM_theta_type))) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" *Malloc error for theta list \n");
		return RADSMM_ERROR_malloc;
	}

	header->theta_count = count;

	memcpy( header->theta_list, list, count*sizeof(RADSMM_theta_type) );

	return 0;
}

                 
int RADSMM_setup_LC( RADSMM_header_type *header, long count )
{

	if (header->model_type != 'D') {
		if ( header->debug_level & 3 )
		    printf(" *Model_type not correct for setup_LC\n ");
		return RADSMM_ERROR_wrong_model;
	}

	if ( (count <= 0) || (count > RADSMM_MAX_liability_classes) ) {
		if ( header->debug_level & 3 )
		    printf(" *Count for setup_LC bad %li\n",count);
		return RADSMM_ERROR_badparam;
	}
	header->LC_count = count;

	return 0;
}

int RADSMM_setup_penetrance( RADSMM_header_type *header, 
			     int    LC_ndx,
                             float  plist1[], 
			     float  plist2[], 
			     float  plist3[], 
			     long    count )
{
/* note that the sent list is in A1,A2,A3...B1,B2,B3... C1,C2,C3 order
   I want to store data in A1,B1,C1,A2,B2,C2,A3,B3,C3,... */

	long	i;
	RADSMM_penetrance_type		*penptr;

	if (header->model_type != 'D') {
		if ( header->debug_level & 3 )
		    printf(" *Model_type not correct for setup_penet\n ");
		return RADSMM_ERROR_wrong_model;
	}

	if ( (count <= 0) || (count>RADSMM_MAX_penetrances) ) {
		if ( header->debug_level & 3 )
		    printf(" *Count for setup_penetrances bad %li\n",count);
		return RADSMM_ERROR_badparam;
	}

	if ( (LC_ndx < 0) || (LC_ndx >= header->LC_count) ) {
		if ( header->debug_level & 3 )
		    printf(" *LC index for setup_penetrances is bad %i\n",LC_ndx);
		return RADSMM_ERROR_badparam;
	}

	if ((header->penetrance_list[LC_ndx] = (RADSMM_penetrance_type *) malloc(count*
                                            sizeof(RADSMM_penetrance_type))) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" *Malloc error for penetrance list \n");
		return RADSMM_ERROR_malloc;
	}

	header->penetrance_count = count;

	penptr = header->penetrance_list[LC_ndx];

	for ( i=0 ; i<count ; i++ ) {
	   penptr->penetrance[0] = plist1[i];
	   penptr->penetrance[1] = plist2[i];
	   penptr->penetrance[2] = plist3[i];
	   penptr++;
	}

	if ( header->debug_level & 16 ) {
		printf( " ***penetrance cnt %li\n", header->penetrance_count );
	}
	return 0;
}



/* Changed double list[] to RADSMM_diseq_type list[] - ams */
int RADSMM_setup_diseq( RADSMM_header_type *header, 
			   RADSMM_diseq_type   list[], 
			   long  count )
{

	if ( header->use_Diseq == 'N' ) {
		if ( header->debug_level & 3 )
		    printf(" *use_diseq flag not correct for setup_diseq\n ");
		return RADSMM_ERROR_wrong_model;
	}

	if ( (count<=0) || (count>RADSMM_MAX_diseqs) ) { 
		if ( header->debug_level & 3 )
		    printf(" *Count for setup_diseq bad %li\n",count);
		return RADSMM_ERROR_badparam;
	}

	if ((header->diseq_list = malloc(count*sizeof(RADSMM_diseq_type))) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" *Malloc error for diseq list \n");
		return RADSMM_ERROR_malloc;
	}

	header->diseq_count = count;

	memcpy( header->diseq_list, list, count*sizeof(RADSMM_diseq_type));

	if ( header->debug_level & 16 ) {
		printf( " ***diseq cnt %li\n", header->diseq_count );
	}
	return 0;
}

int RADSMM_setup_geneFreq( RADSMM_header_type *header, 
			   double   list[], 
			   long  count )
{

	if ( (count<=0) || (count>RADSMM_MAX_geneFreqs) ) { 
		if ( header->debug_level & 3 )
		    printf(" *Count for setup_genefreq bad %li",count);
		return RADSMM_ERROR_badparam;
	}

	if ((header->geneFreq_list = malloc(count*sizeof(RADSMM_geneFreq_type))) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" *Malloc error for genefreq list \n");
		return RADSMM_ERROR_malloc;
	}

	header->geneFreq_count = count;

	memcpy( header->geneFreq_list, list, count*sizeof(RADSMM_geneFreq_type));

	if ( header->debug_level & 16 ) {
		printf( " ***genefreq cnt %li\n", header->geneFreq_count );
	}
	return 0;
}

int RADSMM_setup_qmodel( RADSMM_header_type *header, 
			   RADSMM_quantitative_type   list[], 
			   long  count )
{

	if ( header->model_type != 'Q' ) {
		if ( header->debug_level & 3 )
		    printf(" *Model_type not correct for setup_qmodel\n ");
		return RADSMM_ERROR_wrong_model;
	}

	if ( (count<=0) || (count>RADSMM_MAX_qmodels) ) { 
		if ( header->debug_level & 3 )
		    printf(" *Count for setup_qmodel bad %li",count);
		return RADSMM_ERROR_badparam;
	}

	if ((header->qmodel_list = malloc(count*sizeof(RADSMM_quantitative_type) ) ) == NULL ) {
		if ( header->debug_level & 3 )
		    printf(" *Malloc error for Qmodel list \n");
		return RADSMM_ERROR_malloc;
	}

	header->qmodel_count = count;

	memcpy( header->qmodel_list, list, count*sizeof(RADSMM_quantitative_type) );

	if ( header->debug_level & 16 ) {
		printf( " ***qmodel cnt %li\n", header->qmodel_count );
	}
	return 0;
}


int RADSMM_setup_data( RADSMM_header_type *header, 
				char  dtype,
				int   checksum_type )
{
	if ( header->debug_level & 16 ) {
		printf( " ***entering setup_data \n" );
	}

	if ( (dtype =='f') || (dtype=='F') ) {
		header->element_data_type = 'F';
		header->chunk_size = sizeof( float );
	} else if ( (dtype =='d') || (dtype=='D') ) {
		header->element_data_type = 'D';
		header->chunk_size = sizeof( double );
	} else {
		if (header->debug_level & 3 )
		    printf(" Invalid data type in _setup_data \"%1s\"\n", &dtype );
		return RADSMM_ERROR_badparam;
	}

/*	if ( (header->debug_level & 16) !=0  ) {
		printf( " ***Marker Count is %i %s ", 
			header->marker_count, &header->element_data_type );
		printf( "  Chunksize is %i \n",header->chunk_size);
	}  */

/*	if ( (checksum_type < 0) || (checksum_type > 0) ) {
		if (header->debug_level & 3 )
		    printf(" Invalid checksum_type value %i\n",checksum_type);
		return RADSMM_ERROR_badparam;
	} */
	header->checksum_type = 0;

	return 0;
}

int RADSMM_setup_comments( RADSMM_header_type *header, 
			   char  *comment )
{
#include <time.h>
	time_t  tic;
	char timestring[20];

	if ( comment != NULL )	
	     strncpy( header->description, comment, RADSMM_descr_length );

	if ( header->debug_level & 16 ) {
		printf( " *** comment is '%s'\n", header->description );
	}
	tic = time( &tic );

	strftime( timestring , 18, "%D %H:%M:%S\0", localtime(&tic) );
        /* for some reason, It needs an extra character. */

	strncpy( header->date_string, timestring, 17 );
	
	if ( header->debug_level & 16 ) {
		printf( " *** timestring \"%.17s\"\n", (char *)&header->date_string );
	}
	return 0;
}

int RADSMM_setup_ordering( RADSMM_header_type *header, char order )

{
	switch ( order ) {
	   case 'A':
	   case 'B':
	   case 'C':
	   case 'D': 
	   case 'E': 
	      header->ordering = order;
	      break;
	   default :
		return RADSMM_ERROR_badparam;
	}
	return 0;

}


int RADSMM_file_size( RADSMM_header_type *header, 
		      double *size, 
                      int *number_of_files )
{
#include <limits.h>
#include "RADSMM_file_header.h"
	double  dsize;
	double	datasize;
	long 	fcnt;

/*	header->debug_level = header->debug_level & 255; */

	switch ( header->model_type ) {
	    case 'D':
		dsize = header->pedigree_count 
		      * header->penetrance_count
		      * header->geneFreq_count
		      * header->marker_count;
		break;
	    case 'Q':
		dsize = header->pedigree_count 
		      * header->qmodel_count
		      * header->marker_count;
		break;
	    default:
		if (header->debug_level & 3 )
		    printf(" Invalid model type  \"%1s\"\n", &header->model_type );
		return RADSMM_ERROR_wrong_model;
	}

        if ( header->marker_type == '2' ) {
	   if ( header->theta_matrix_type=='G' )
		dsize = dsize * header->theta_count * header->theta_count;
	   else if ( header->theta_matrix_type=='D' )
		dsize = dsize * header->theta_count;
	   else {
		if (header->debug_level & 3 )
		    printf(" Invalid theta matrix type _file_size \"%1s\"\n", &header->theta_matrix_type );
		return RADSMM_ERROR_badparam;
	   }
        }

        if (header->use_Diseq != 'N' ) {
	   dsize = dsize * header->diseq_count;
        }

	/* now add in headersize */

	datasize = dsize;

	dsize = dsize * header->chunk_size;
	dsize = dsize + sizeof(RADSMM_file_header_type); 
	dsize = dsize + header->pedigree_count * sizeof(RADSMM_pedigree_type);
	dsize = dsize + header->marker_count * header->markerlabel_size;
	dsize = dsize + header->pedigree_count * header->pedigreelabel_size;

	switch ( header->marker_type ) {
	    case '2':
	       dsize = dsize + header->marker_count * sizeof(RADSMM_marker_type);
	       dsize = dsize + header->theta_count * sizeof(RADSMM_theta_type);
               break;
	    case 'M':
	       dsize = dsize + header->marker_count * sizeof(RADSMM_marker_type);
               break;
        }

	switch ( header->model_type ) {
	    case 'D':
		dsize = dsize + header->penetrance_count * sizeof(RADSMM_penetrance_type);
  	        dsize = dsize + header->geneFreq_count * sizeof(RADSMM_geneFreq_type);
		break;
	    case 'Q':
		dsize = dsize + header->qmodel_count * sizeof(RADSMM_quantitative_type);
		break;
	    default:
		if (header->debug_level & 3 )
		    printf(" Invalid model type  \"%1s\"\n", &header->model_type );
		return RADSMM_ERROR_wrong_model;
	}

	switch ( header->use_Diseq ) {
	    case 'N':
               break;
            default :
               dsize = dsize + header->diseq_count * sizeof(RADSMM_diseq_type);
        }



	if ( header->debug_level & 16)
		printf("*** file size is %f\n",dsize);

	if ( dsize >= LONG_MAX ) {
		/* here is where we determine # of files and file size */
  		fcnt = 1 + (datasize*header->chunk_size / (double) LONG_MAX); /* truncation occurs */	
		header->chunks_per_file = datasize / (double)fcnt;
		if ( fcnt * header->chunks_per_file < datasize ) header->chunk_size++;
		*size = header->chunks_per_file * header->chunk_size;
		header->number_of_files = fcnt;
		*number_of_files = header->number_of_files+1;
		if ( header->debug_level & 16 )
		    printf( " ***Multi files %i * %li(%E) = %E\n", header->number_of_files,
					    header->chunks_per_file,*size,datasize );
		return 0;
	} else if ( dsize <= 0 ) {
		if ( header->debug_level & 3 ) 
		    printf("   for some reason the file size is <=0, %f",dsize);
		*size = -1;
		return RADSMM_ERROR_badparam;
	} else {
		*size = dsize;
		header->chunks_per_file = datasize;
		header->number_of_files = 0;
		*number_of_files = 1;
	}

	return 0;
}
