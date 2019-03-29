#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "RADSMM.h"

#define TST_MARKTYPE '2'
#define TST_MODEL    'Q'
#define TST_USEDISEQ 'N'
#define TST_MARKERS  12
#define TST_PEDS     30
#define TST_THETA    30
#define TST_PENS     1 /*36*/
#define TST_GENES    1 /*10*/
#define TST_DISEQ    1 /*10*/
#define TST_QMOD     5
#define THETA_TYPE   'D'


int  fakeval( long pedval, long thetaval, long penval, long geneval, 
	      long marker, double *buff )
/* Calculate fake values for data to write, We'll use the index parameters to get
unique values. */
{

	*buff =  (pedval) + (thetaval*100.0) + (penval*1000.0) +
		      (geneval*10000.0) + (marker*0.01);
	return 0;
}


long ranrange( long low, long high )
/* generate a random integer between the high and low values. */
{
	float f;
	f = (float) random() / (float)RAND_MAX;
	f = low + f * (high-low+1);
	if ( f >= high+1.0 ) f = high;
	return (long) f;
}

int ranseq( int array[], int size )
/* fill array with random sequence of 0 to size-1; */
{
	int i, j, val;

	for (i=0; i <size ; i++ ) {
		array[i] = i;
	}
	
	for( i=0 ; i<size-1 ; i++ ) {
		j = ranrange( i, size-1 );
		val = array[j];
		array[j] = array[i];
		array[i] = val;
	}
/*	for (i=0;i<size;i++) printf( " %i,",array[i]);
	printf("\n");
*/	return 0;
}


int main( int argc, char *argv[] )
{


        RADSMM_marker_type      marks[TST_MARKERS] = 
                 { 0.0, 1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 
                   6.0, 7.0, 8.0, 9.0  }; 
/*	RADSMM_pedigree_type 	peds[TST_PEDS] =
		{ 1,2,3,4,5, 10,11,12,13,14, 20,21,22,23,24, 30,31,32,33,34,  
		  41,44,46,47,49 }; */
	RADSMM_theta_type	thetas[TST_THETA] =
		{ .01, .02, .03, .04, .05, .06, .07, .08, .09, .10,
		  .11, .12, .13, .14, .15, .16, .17, .18, .19, .20,
		  .21, .22, .23, .24, .25, .26, .27, .28, .29, .30 };
/*		  .31, .32, .33, .34, .35, .36, .37, .38, .39, .40,
		  .41, .42, .43, .44, .45, .46, .47, .48, .49, .50 }; */

        static float            pens1[TST_PENS] =
         { 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1, 0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2
          ,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4, 0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6 };
        static float            pens2[TST_PENS] =
         { 0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9,  0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9,
           0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9,  0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9 };
        static float            pens3[TST_PENS] =
         { 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7,
           0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7 };

        RADSMM_quantitative_type qmods[TST_QMOD] = {
		{0.1, 0.1, 0.2, 0.1, 0.2, 0.1},
		{0.1, 0.2, 0.2, 0.2, 0.2, 0.1},
		{0.1, 0.2, 0.3, 0.4, 0.2, 0.1},
		{0.1, 0.3, 0.4, 0.1, 0.2, 0.3},
		{0.2, 0.1, 0.1, 0.1, 0.3, 0.1} };


	RADSMM_geneFreq_type 	genes[TST_GENES] = 
	  { .010,.050,.100,.125,.150,.175,.200,.225,.250,.275 };
/* 	    .300,.310,.320,.330,.340,.350,.375,.400,.425,.450 }; */

	char	comment[29] = "Test data file gened by tst2";

	/* regular vars */	

	long 	pedndx, thetandx, penndx, gfndx, mrkndx, qndx, disndx;
	long    i, i2;
	long long j;
	int	stat;
	double 	fsize;
	int	fcnt;
	RADSMM_header_type hdr;
	RADSMM_header_type hdr2;
	char	filename[] = "testdatafile.dat";
	char	label[23];
	double	data, test;
	double	cputime, icputime;
	time_t	walltime, iwalltime; 
/*	int	order1[999];
	int	order2[999];
*/
	printf(" Starting...\n\n");

	stat = RADSMM_setup_init( &hdr , 127 );
	if ( stat != 0 ) {
	    printf(" Error: setup_init %i\n", stat );
	    exit(2);
	}

        stat = RADSMM_setup_type( &hdr, TST_MARKTYPE, TST_MODEL, TST_USEDISEQ );
        if ( stat != 0 ) {
            printf(" Error: setup type %i\n",stat );
            exit(2);
        }

	stat = RADSMM_setup_ordering( &hdr, 'A' );
	if ( stat != 0 ) {
	    printf(" Error: setup ordering %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_marker( &hdr, marks ,(long) TST_MARKERS );
	if ( stat != 0 ) {
	    printf(" Error: setup peds %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_pedigree( &hdr, NULL ,(long) TST_PEDS );
	if ( stat != 0 ) {
	    printf(" Error: setup peds %i\n",stat );
	    exit(2);
	}

	if ( TST_MARKTYPE == '2') {
           stat = RADSMM_setup_theta( &hdr, thetas, (long) TST_THETA, THETA_TYPE );
           if ( stat != 0 ) {
	      printf(" Error: setup theta %i\n",stat );
	      exit(2);
           }
	}

        if ( TST_MODEL == 'D' ) {
           stat = RADSMM_setup_LC( &hdr, 3 );
           if ( stat != 0 ) {
              printf(" Error: setup_LC %i\n",stat );
              exit(2);
           }

	   stat = RADSMM_setup_penetrance( &hdr, 0, pens1,pens2,pens3, (long) TST_PENS );
	   if ( stat != 0 ) {
	      printf(" Error: setup penet %i\n",stat );
	      exit(2);
	   }
	   stat = RADSMM_setup_penetrance( &hdr, 1, pens2,pens3,pens1, (long) TST_PENS );
	   if ( stat != 0 ) {
	      printf(" Error: setup penet %i\n",stat );
	      exit(2);
	   }
	   stat = RADSMM_setup_penetrance( &hdr, 2, pens3,pens1,pens2, (long) TST_PENS );
	   if ( stat != 0 ) {
	      printf(" Error: setup penet %i\n",stat );
	      exit(2);
	   }

	   stat = RADSMM_setup_geneFreq( &hdr, genes, (long) TST_GENES );
	   if ( stat != 0 ) {
	     printf(" Error: setup genefreq %i\n",stat );
	     exit(2);
	  }

        } else {
           
	   stat = RADSMM_setup_qmodel( &hdr, qmods, (long) TST_QMOD );
	   if ( stat != 0 ) {
	      printf(" Error: setup qmodels %i\n",stat );
	      exit(2);
	   }
        }
      
        if ( TST_USEDISEQ == 'Y' ) {
	   stat = RADSMM_setup_diseq( &hdr, genes, (long) TST_GENES );
	   if ( stat != 0 ) {
	      printf(" Error: setup diseq %i\n",stat );
	      exit(2);
	   }
        }

	stat = RADSMM_setup_data( &hdr,'D',0 );
	if ( stat != 0 ) {
	    printf(" Error: setup data %i\n",stat );
	    exit(2);
	}
	printf( " ---data type %1s---\n", &hdr.element_data_type );

	stat = RADSMM_setup_comments( &hdr, comment );
	if ( stat != 0 ) {
	    printf(" Error: setup comments %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_markerlabel( &hdr, (int)12 );
	if ( stat != 0 ) {
	    printf(" Error: setup comments %i\n",stat );
	    exit(2);
	}

	for ( i=1 ; i <= TST_MARKERS ; i++ ) {
		sprintf( label, " MRK - %li",i );
		stat = RADSMM_set_markerlabel( &hdr, i-1, label );
		if ( stat != 0 ) {
		    printf(" Error: set markerlabel %i\n",stat );
		    exit(2);
		}
	}

	stat = RADSMM_setup_pedigreelabel( &hdr, (int)20 );
	if ( stat != 0 ) {
	    printf(" Error: setup comments %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_file_size( &hdr, &fsize, &fcnt );
	if ( stat != 0 ) {
	    printf(" Error: file size (%f)%i\n",fsize, fcnt );
	    exit(2);
	}
	
	printf(" File size is %f * %i\n",fsize, fcnt );
	system( "date" );

	stat = RADSMM_create_file( &hdr, filename, 0777 );
	if ( stat != 0 ) {
	    printf(" Error: createfile %i\n", stat );
	    exit(2);
	}

	printf( " ---data type %c---\n", hdr.element_data_type );
	printf(" File Created. \"%c\" \n",hdr.open_flag );


/* now fill the whole thing */

	j = 0;
	icputime = cputime = (double)clock() / CLOCKS_PER_SEC;
	iwalltime = walltime = time(NULL);

	system( "date" );

/*	ranseq( order2, TST_GENES ); */
/*  A ordering: Pedigree, Marker, MaleTheta, FemaleTheta, GeneFreq, Penetrance, Qmodel, diseq */


	for ( pedndx=0 ; pedndx<TST_PEDS ; pedndx++ ) {
           for ( mrkndx=0 ; mrkndx<TST_MARKERS ; mrkndx++ ) {
              if (THETA_TYPE == 'G' ) 
                 i2 = (TST_THETA * TST_THETA);
              else
                 i2 = TST_THETA;
              for ( thetandx=0 ; thetandx<i2 ; thetandx++ ) {
	         for ( gfndx=0 ; gfndx<TST_GENES ; gfndx++ ) {
                    for ( penndx=0 ; penndx<TST_PENS ; penndx++ ) {
                       for ( qndx=0 ; qndx<TST_QMOD ; qndx++ ) {
                          for ( disndx=0 ; disndx<TST_DISEQ ; disndx++ ) {
	
 				stat = fakeval( pedndx, thetandx, penndx, 
						gfndx, mrkndx, &data );

	/*	printf(" data = %E\n",data ); */

                                stat = RADSMM_write_data( &hdr, pedndx, mrkndx, 
                                   thetandx, gfndx, penndx, qndx, disndx , &data );
				if ( stat != RADSMM_ERROR_success ) {
	    	  		  printf(" Error: writing (#%lli): err %i\n",j,stat );
				  printf( " debug=%i, open=%c\n",hdr.debug_level, hdr.open_flag);
	    	  		  exit(2);
				}

		/* read the data back and see if it's different */
		/*
                                stat = RADSMM_read_data( &hdr2, pedndx, mrkndx, 
                                     thetandx, gfndx, penndx, qndx, disndx, &test );
				if ( stat != 0 ) {
	    			  printf(" Error: read (#%i): err %i\n",j,stat );
	    			  exit(2);
				}

				if (test == RADSMM_EMPTY ) {
			 	  printf ( " %i hmm empty %E\n",j,data );
				} else if ( abs(data-test) > data*0.000001 ) { 
		  		  printf ( "  bad echo =  %i,%i,%i,%i,%i  %E(%E) %E\n",
		 		  	pedndx,thetandx,penndx,gfndx,mrkndx,
					data,test,(test-data));
				} */
				j++;
				if ( (j % 20000L) == 0 )
				    printf("*");
		                if ( (j % 1000000L) == 0 ) {
        	                   printf(" %f, %f\n",((double)clock()/CLOCKS_PER_SEC)
					-cputime,difftime(time(NULL),walltime) );
				   walltime = time(NULL);
				   cputime = clock()/CLOCKS_PER_SEC;

				}			    
                          }
                       }
                    }
                 }
	      }
	   }
	}
	printf("\n");

	system( "date" );
	printf(" %lli points\n",j);
	printf(" Cputime per 1M= %f\n", (((double)clock()/CLOCKS_PER_SEC)-icputime)*1000000.0/j );
	printf(" Walltime per 1M= %f\n", difftime(time(NULL),iwalltime)*1000000.0/j );


	stat = RADSMM_close_file( &hdr );
	
	stat = RADSMM_setup_init( &hdr2, 255 );
	stat = RADSMM_open_file( &hdr2, filename, 'r' );
	if ( stat != 0 ) {
		printf( " ERROR opening file: %i\n",stat );
		exit(2);
	}
	for ( i = 0 ; i < hdr2.marker_count ; i++ )
	    printf( " Marker %li label \"%s\"\n",i,RADSMM_get_markerlabel(&hdr2,i));

	for ( i = 0 ; i < hdr2.pedigree_count ; i++ )
	    printf( " Pedigree %li label \"%s\"\n",i,RADSMM_get_pedigreelabel(&hdr2,i));


/* read them back randomly and check */

	printf("\n Reading back and checking \n");
	srandom( 12345 ); 
	for (j=0;j<8000000;j++){
		pedndx = ranrange(1,TST_PEDS-1);
		if ( pedndx < 0 ) {
		    printf(" Error: index ped (%lli) %li\n",j,pedndx );
		    exit(2);
		}

		if (THETA_TYPE == 'G' ) 
			thetandx  = ranrange(0,(TST_THETA-1)*(TST_THETA-1) );
		else
			thetandx  = ranrange(0,(TST_THETA-1) );

		penndx = ranrange(0,TST_PENS-1);

		gfndx = ranrange(0,TST_GENES-1);

		mrkndx = ranrange( 0, TST_MARKERS-1 );

                qndx = ranrange(0, TST_QMOD-1);
                disndx = ranrange(0, TST_DISEQ-1);

		stat = fakeval( pedndx, thetandx, penndx, gfndx, mrkndx, 
						&test );
                stat = RADSMM_read_data( &hdr2, pedndx, mrkndx, 
                                     thetandx, gfndx, penndx, qndx, disndx, &data );

		if ( stat != 0 ) {
	    	  printf(" Error: read (#%lli): err %i\n",j,stat );
	    	  exit(2);
		}

		if ( data == RADSMM_EMPTY ) {
			printf ( "%lli hmm empty %E\n",j,test );
		} else if ( abs(data-test) > test*.000001 ) { 
		  	printf ( "%lli not equal = %li,%li,%li,%li,%li  %E(%E)\n",
			  j,pedndx,thetandx,penndx,gfndx,mrkndx,data,test );
		}

		if ( (j % 7000L) == 0 ) 
			putchar( (int)'*' );

		if ( (j % 490000L) == 0 ) printf(" %lli\n",j);

	}

	printf("\n");

	stat = RADSMM_close_file( &hdr2 );

	return 0;

}
