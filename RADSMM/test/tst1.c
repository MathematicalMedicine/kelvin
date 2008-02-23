#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../RADSMM.h"

#define TST_MARKERS  40
#define TST_MARKTYPE '2'
#define TST_PEDS     50
#define TST_MODEL    'D'
#define TST_THETA    50
#define TST_PENS     36
#define TST_GENES    20
#define THETA_TYPE 'D'
#define TST_DISEQ 'N'


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
	f = f * (high-low+1);
	return (long) f + low;
}


int main( int argc, char *argv[] )
{

/*	pedigree_type 	peds[TST_PEDS] =
		{ 1,2,3,4,5, 6,7,8,9,10,
		11,12,13,14,15,16,17,18,19,20,
		21,22,23,24,25,26,27,28,29,30,
		31,32,33,34,35,36,37,38,39,40,
		41,42,43,44,45,46,47,48,49,50,
		51,52,53,54,55,56,57,58,59,60 }; */
	RADSMM_theta_type	thetas[TST_THETA] =
		{ .01, .02, .03, .04, .05, .06, .07, .08, .09, .10,
		  .11, .12, .13, .14, .15, .16, .17, .18, .19, .20,
		  .21, .22, .23, .24, .25, .26, .27, .28, .29, .30,
		  .31, .32, .33, .34, .35, .36, .37, .38, .39, .40,
		  .41, .42, .43, .44, .45, .46, .47, .48, .49, .50 };
        static float            pens1[TST_PENS] =
         { 0.1,0.1,0.1 ,0.1,0.1,0.1, 0.1,0.1,0.1, 0.2,0.2,0.2, 0.2,0.2,0.2, 0.2,0.2,0.2
          ,0.4,0.4,0.4, 0.4,0.4,0.4, 0.4,0.4,0.4, 0.6,0.6,0.6, 0.6,0.6,0.6, 0.6,0.6,0.6 };
        static float            pens2[TST_PENS] =
         { 0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9,  0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9,
           0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9,  0.1,0.1,0.1, 0.3,0.3,0.3, 0.9,0.9,0.9 };
        static float            pens3[TST_PENS] =
         { 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7,
           0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7, 0.0,0.5,0.7 };

	RADSMM_geneFreq_type 	genes[TST_GENES] = 
	  { .010,.050,.100,.125,.150,.175,.200,.225,.250,.275,
 	    .300,.310,.320,.330,.340,.350,.375,.400,.425,.475 };

	char	comment[17] = "Test data file";

	/* regular vars */	

	long 	pedndx, thetandx, penndx, gfndx, mrkndx;
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
	double 	cputime, icputime;
	double 	val1, val2, val3;
	time_t	walltime, iwalltime;

	printf(" Starting...\n\n");

	stat = RADSMM_setup_init( &hdr , 255 );
	if ( stat != 0 ) {
	    printf(" Error: setup_init %i\n", stat );
	    exit(2);
	}

        stat = RADSMM_setup_type( &hdr, TST_MARKTYPE, TST_MODEL, TST_DISEQ );
	if ( stat != 0 ) {
	    printf(" Error: setup type %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_ordering( &hdr, 'B' );
	if ( stat != 0 ) {
	    printf(" Error: setup ordering %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_pedigree( &hdr, NULL ,(long) TST_PEDS );
	if ( stat != 0 ) {
	    printf(" Error: setup peds %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_theta( &hdr, thetas, (long) TST_THETA, THETA_TYPE );
	if ( stat != 0 ) {
	    printf(" Error: setup theta %i\n",stat );
	    exit(2);
	}

        stat = RADSMM_setup_LC( &hdr, 1 );
	if ( stat != 0 ) {
	    printf(" Error: setup_LC %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_penetrance( &hdr, 0, pens1,pens2,pens3, (long) TST_PENS );
	if ( stat != 0 ) {
	    printf(" Error: setup penet %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_geneFreq( &hdr, genes, (long) TST_GENES );
	if ( stat != 0 ) {
	    printf(" Error: setup genefreq %i\n",stat );
	    exit(2);
	}

	stat = RADSMM_setup_data( &hdr, (int)TST_MARKERS,'F',0 );
	if ( stat != 0 ) {
	    printf(" Error: setup data %i\n",stat );
	    exit(2);
	}
	printf( " ---data type %s---\n", &hdr.element_data_type );

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
		sprintf( label, " MRK - %li\0",i );
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
	    printf(" Error: file size (%f)%i  %i\n",fsize, fcnt, stat );
	    exit(2);
	}
	
	printf(" File size is %f * %i\n",fsize, fcnt );

	stat = RADSMM_create_file( &hdr, filename, 0777 );
	if ( stat != 0 ) {
	    printf(" Error: createfile %i\n", stat );
	    exit(2);
	}

	printf( " ---data type %s---\n", &hdr.element_data_type );
	printf(" File Created.\n");

/* now fill randomly */
	srandom( 12345 ); 

	icputime = cputime = (double)clock() / CLOCKS_PER_SEC;
	iwalltime = walltime = time(NULL);
	for (j=0;j<=1000000;){
		pedndx = ranrange(1, TST_PEDS-1);


		if (THETA_TYPE == 'G' ) 
			val2 = TST_THETA*TST_THETA;
		else
			val2 = TST_THETA;
		thetandx = ranrange(0,val2-1);
/*		val1 = val1 + (( ((float)random()/(float)RAND_MAX)-0.5)*2.2E-5);
		val2 = val2 + (( ((float)random()/(float)RAND_MAX)-0.5) *val2*1.0E-5);
		thetandx = RADSMM_index_theta( &hdr, val1, val2 );
		if ( thetandx < 0 ) {
		    printf(" Error: index theta (%i) %i,%i : %i \n",j,i,i2,thetandx );
		    exit(2);
		} */

		i = ranrange(0,TST_PENS-1);
		val1 = pens1[i];
		val2 = pens2[i];
		val3 = pens3[i];
		penndx = RADSMM_index_penetrance( &hdr, 1, (float)val1, (float)val2, (float)val3 );
		if ( penndx < 0 ) {
		    printf(" Error: index penet (%lli) %li,%li\n",j,i,penndx );
		    exit(2);
		}

		i = ranrange(0,TST_GENES-1);
		val1 = genes[i];
		gfndx = RADSMM_index_genefreq( &hdr, val1 );
		if ( gfndx < 0 ) {
		    printf(" Error: index gene freq (%lli) %li, %li\n",j,i,gfndx );
		    exit(2);
		} 

		mrkndx = ranrange( 0, TST_MARKERS-1 );

		stat = fakeval( pedndx, thetandx, penndx, gfndx, mrkndx, 
								&data );

	/*	printf(" data = %E\n",data ); */

		stat = RADSMM_write_data( &hdr, pedndx, mrkndx, thetandx, gfndx, penndx,
						1, 1 , &data );
		if ( (stat != 0) && ( stat != RADSMM_ERROR_writeover_valid_data) ) {
	    	  printf(" Error: writing (#%lli): err %i, %e  \n",j,stat,data );
	    	  exit(2);
		}

		if ( (j % 50000L) == 0 ) {
			printf("%lli %f, %f\n",j,((double)clock()/CLOCKS_PER_SEC)-cputime,difftime(time(NULL),walltime) );
			walltime = time(NULL);
			cputime = clock()/CLOCKS_PER_SEC;
		}
		j++;
	}

	printf("\n");

        printf(" Cputime per 1M= %f\n", (((double)clock()/CLOCKS_PER_SEC)-icputime)*1000000.0/j );
        printf(" Walltime per 1M= %f\n", difftime(time(NULL),iwalltime)*1000000.0/j );


	stat = RADSMM_close_file( &hdr );
	

	stat = RADSMM_setup_init( &hdr2, 255 );
	stat = RADSMM_open_file( &hdr2, filename, 'r' );
	if ( stat != 0 ) {
		printf( " ERROR opening file: %i\n",stat );
		exit(2);
	}

/* now, reset the random sequence and read them back and check */

	j = 0;
	printf("\n Reading back and checking \n");

        for ( pedndx=0 ; pedndx<TST_PEDS ; pedndx++ ) {
                                  
          for ( mrkndx=0 ; mrkndx<TST_MARKERS ; mrkndx++ ) {
                                   
            if (THETA_TYPE == 'G' )
              i2 = (TST_THETA * TST_THETA);
            else
              i2 = TST_THETA;
            for ( thetandx=0 ; thetandx<i2 ; thetandx++ ) {
                                  
              for ( gfndx=0 ; gfndx<TST_GENES ; gfndx++ ) {
        
                for ( penndx=0 ; penndx<TST_PENS ; penndx++ ) {
                                 
		stat = fakeval( pedndx, thetandx, penndx, gfndx, mrkndx, 
						&test );
		stat = RADSMM_read_data( &hdr2, pedndx, mrkndx, thetandx, gfndx, penndx,
						1, 1, &data );
		if ( stat != 0 ) {
	    	  printf(" Error: read (#%lli): err %i\n",j,stat );
	    	  exit(2);
		}

		if ( data == RADSMM_EMPTY ) {
			/* printf ( "%i hmm empty %E\n",j,test ); */
		} else if ( abs(data-test) > test*.000001 ) { 
		  	printf ( "%lli not equal = %li,%li,%li,%li,%li  %E(%E)\n",
			  j,pedndx,thetandx,penndx,gfndx,mrkndx,data,test );
		}
		j++;
		if ( (j % 100000L) == 0 ) 
			putchar( (int)'*' );
		/*	printf("*"); */

		if ( (j % 1000000L) == 0 ) printf(" %lli\n",j);
		}
	      }
	    }
	  }
	}

	printf("\n");

	stat = RADSMM_close_file( &hdr2 );

	return 0;

}
