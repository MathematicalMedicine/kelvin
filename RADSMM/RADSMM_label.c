/* This contains some of the marker & pedigree labeling  routines */
#include <stdio.h>
#include <errno.h>
#include <malloc.h>
#include "RADSMM.h"

int RADSMM_setup_markerlabel (	RADSMM_header_type *header, int label_length ) 
/* check the parameter for sanity, allocate memory */
{
	int	i;
	char	*p;
	size_t 	mlen;


	if ( label_length == 0 ) {
	    	header->markerlabel_size = label_length;
		header->markerlabel_list = NULL;
		return 0;
        } else if ( (label_length<0) || (label_length>RADSMM_MAX_markerlabel_length) ) {
                if ( header->debug_level & 3 )
                    printf(" Count for setup_markerlabel bad %i",label_length);
                return RADSMM_ERROR_badparam;
        }

        if ((header->markerlabel_list = malloc(label_length*header->marker_count*sizeof(char) ) ) == NULL ) {
                if ( header->debug_level & 3 )
                    printf(" *Malloc error for pedigree list. %i\n",errno);
                return RADSMM_ERROR_malloc;
        }

	if ( header->debug_level &16 )
	   printf(" *** Markerlabel space allocated %i,%li\n",label_length,header->marker_count);

	/* Initialize the markerlabel depeing on size */
        header->markerlabel_size = label_length;
	mlen = header->markerlabel_size;
	p = header->markerlabel_list;	
	for ( i=1 ; i <= header->marker_count ; i++ ) {
		if ( header->markerlabel_size == 1  )
			snprintf( p, mlen, "%1i", i & 10 );
		else if ( header->markerlabel_size == 2 ) 
			snprintf( p, mlen, "%2i", i & 100 );
		else if ( header->markerlabel_size == 3 )
			snprintf( p, mlen, "%3i", i  );
		else if ( header->markerlabel_size <= 8 )
			snprintf( p, mlen, "M%i", i  );
		else 
			snprintf( p, mlen, "Marker%i", i );
		p = p+header->markerlabel_size;
	}


        return 0;
}


int RADSMM_set_markerlabel( RADSMM_header_type *header, 
			    int markerindex,
                            char *label )
/* set the label in memory */
{

	char *src;
	char *dest;
	int	i;


	if ( (markerindex >= header->marker_count) || (markerindex < 0) ) {
		if ( header->debug_level & 3 )
		     printf(" Marker Label index out of range %i \n",markerindex);
		return RADSMM_ERROR_badparam;
	}

	dest = header->markerlabel_list + (header->markerlabel_size*markerindex);
	src = label;

	for ( i=0 ; (i<header->markerlabel_size) && (src!="\0") ; i++)
		*dest++ = *src++;

	return 0;
}


int RADSMM_setup_pedigreelabel ( RADSMM_header_type *header, int label_length ) 
/* check the parameter for sanity, allocate memory */
{
	int	i;
	char	*p;
	size_t	plen;

	if ( label_length == 0 ) {
	    	header->pedigreelabel_size = label_length;
		header->pedigreelabel_list = NULL;
		return 0;
        } else if ( (label_length<0) || (label_length>RADSMM_MAX_pedigreelabel_length) ) {
                if ( header->debug_level & 3 )
                    printf(" Count for setup_pedigreelabel bad %i",label_length);
                return RADSMM_ERROR_badparam;
        }

        if ((header->pedigreelabel_list = malloc(label_length*header->pedigree_count*sizeof(char) ) ) == NULL ) {
                if ( header->debug_level & 3 )
                    printf(" *Malloc error for pedigree list. %i\n",errno);
                return RADSMM_ERROR_malloc;
        }

	/* Initialize the pedigreelabel depeing on size */
        header->pedigreelabel_size = label_length;
	plen = header->pedigreelabel_size;
	p = header->pedigreelabel_list;	
	for ( i=1 ; i <= header->pedigree_count ; i++ ) {
		if ( header->pedigreelabel_size == 1  )
			snprintf( p, plen, "%1i", (i & 10) );
		else if ( header->pedigreelabel_size == 2 ) 
			snprintf( p, plen, "%2i", i & 100 );
		else if ( header->pedigreelabel_size == 3 )
			snprintf( p, plen, "%3i", i  );
		else if ( header->pedigreelabel_size <= 10 )
			snprintf( p, plen, "P%i", i  );
		else 
			snprintf( p, plen, "Pedigree%i", i );
		p = p + header->pedigreelabel_size;
	}


        return 0;
}


int RADSMM_set_pedigreelabel( RADSMM_header_type *header, 
			    int pedigreeindex,
                            char *label )
/* set the label in memory */
{

	char *src;
	char *dest;
	int	i;


	if ( (pedigreeindex >= header->pedigree_count) || (pedigreeindex < 0) ) {
		if ( header->debug_level & 3 )
		     printf(" pedigree Label index out of range %i \n",pedigreeindex);
		return RADSMM_ERROR_badparam;
	}

	dest = header->pedigreelabel_list + (header->pedigreelabel_size*pedigreeindex);
	src = label;

	for ( i=0 ; (i<header->pedigreelabel_size) && (src!="/0") ; i++)
		*dest++ = *src++;

	return 0;
}

char  *RADSMM_get_markerlabel( RADSMM_header_type *header, int marker_index )
/* fetch the markerlabel and return it */

{
	int i;
	char *p, *p2;

	if ( (marker_index < 0) || (marker_index>=header->marker_count) ) {
		if ( header->debug_level & 3 )
		   printf( " IN _get_markerlabel, marker_index out fo range %i\n",marker_index);
		   return "Error\0";
	}

	p = header->markerlabel_list + marker_index*header->markerlabel_size;

	/* check for a null */
	
	for ( i=0,p2=p ; i<header->markerlabel_size ; i++,p2++ )
		if ( *p2 == '\0' ) {
			return p;
		}

	/* No null... Must allocate a string, copy buffer, append null */
	if  (header->debug_level & 16) {
		printf( " ***Allocating for markerlabel %i\n",marker_index);
		printf( " ***\"%20s\"\n", p);
	}

	p2 = calloc( header->markerlabel_size+1, sizeof(char) );
				/* pre initialized to null */
	if ( p2 == NULL ) {
		if (header->debug_level &3 ) 
		  printf(" RADSMM_get_markerlabel:Error allocating for string\n");
		return "Error\0";
	}

	for ( i=0 ; i<header->markerlabel_size ; i++ )
		*(p2+i) = *(p+i);

	return p2;
}

char  *RADSMM_get_pedigreelabel( RADSMM_header_type *header, int ped_index )
/* fetch the markerlabel and return it */

{
	int i;
	char *p, *p2;;

	if ( (ped_index < 0) || (ped_index>=header->pedigree_count) ) {
		if ( header->debug_level & 3 )
		   printf( " IN _get_pedigreelabel, pedigree_index out of range %i\n",ped_index);
		   return "Error\0";
	}

	p = header->pedigreelabel_list + ped_index*header->pedigreelabel_size;

	/* check for a null */
	
	for ( i=0,p2=p ; i<header->pedigreelabel_size ; i++,p2++ )
		if ( *p2 == '\0' ) {
			return p;
		}

	/* No null... Must allocate a string, copy buffer, append null */

	p2 = calloc( header->pedigreelabel_size+1, sizeof(char) );
				/* pre initialized to null */
	if ( p2 == NULL ) {
		if (header->debug_level &3 ) 
		  printf(" RADSMM_get_markerlabel:Error allocating for string\n");
		return "Error\0";
	}

	for ( i=0 ; i<header->pedigreelabel_size ; i++ )
		*p2++ = *p++;

	return p2;
}

