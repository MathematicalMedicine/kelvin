/****************** NflankPick.c ****************/
/* This program is part of the nflank suite. The wrapper program 
read the input data files and send part of that information to
this program where parameters are generated to control the rest of the 
run */

/* Input File Format 
   Line 1  Method #[M] [C#]
	2  step size, Left offset, right offset
	3  K or H for Kosambi or Haldane
	3+ list of markers spacings.
	Last line -1 or EOF
   Output File Format

*/

/* changelog ...
 3-29-05 m&m changed > to >= fixing bug when TrialLoci on top of marker at the
             end of list
*/


#include <stddef.h>
#include <stdio.h>
#include <math.h>

#define MAXLOC 5000
#define GO_LEFT -1
#define GO_RIGHT 1
#define FREE 0
#define USED 1


typedef struct method {
	int EvenCnt;
	int MidFlag;
	int CloseCnt;
	int TotalCnt;
} method_type;


/* Globals */
	double	locations[ MAXLOC ]; /* location of markers */
	int	MarkerCnt;  /* how many entries in locations and in_use */
	int	in_use[ MAXLOC ]; /* set to either FREE or USED */
	method_type	method;
	double	StartOffset, StopOffset, inc;
	char	mapping;  /* either H or K */



double to_theta( double in )   
/* Convert the input centimorgan to a theta value. */
/*   Uses external variable "mapping" to determine if Haldane or Kosambi */
{
   if ( mapping == 'H' ) { /* do Haldane */
      return 0.5*(1.0 - exp(-2.0*fabs(in/100.0) ) );
   } else { /* do Kosambi */
      return  0.5*((exp(4.0*(in/100.0)) - 1.0)/(exp(4.0*(in/100.0)) + 1.0));
   }
}


int where_at( double placeme )
/* Return the index in location the current loc belongs such that
 location[i-1] <= placeme && location[i]>placeme */
{
	int i;

	if ( placeme < 0.0 ) {
	    return -1;  /* should never happen */
	}
	if ( placeme > locations[MarkerCnt] ) {
	    return MarkerCnt+1;  /* should never happen */
	}

	if ( placeme <= locations[0]) {
	    return 0;
	}

	if ( placeme >= locations[MarkerCnt-1]) {
	    return MarkerCnt; /* oddball, Allow last point to work */
	}

	for ( i = 1 ; i <= MarkerCnt ; i++ ) {
	    if ( (placeme >= locations[i-1]) && (placeme < locations[i]) )
		return i;
	}
	fprintf( stderr," Someting messed up in 'where_at', %f\n",placeme );
	fprintf( stderr," Aborting\n");
	exit (-1);

}


int First_Free( int direction , double spot )
/* Return the index to the first free marker Using the in_use array */
/* uses in_use, */
{
	int 	bounce, i, step, start;

	bounce = 0;

	  /* printf( "FF start %3.1f dir %i", spot, direction); */

	start = where_at( spot );

	if ( start == 0 ) 
		direction = GO_RIGHT;
	if ( start == MarkerCnt  )
		direction = GO_LEFT;

	if ( direction == GO_RIGHT ) {
	    step = 1;
	} else { /* GO_LEFT */
	    step = -1;
	    start = start - 1;
	}

	i = start;

	while ( in_use[i] != FREE ) {
	
	    i = i + step;
	    if ( (i < 0) || (i > (MarkerCnt-1) ) ) { /* arrg, at the start or end*/
		bounce++;
		if ( bounce >= 2) {
		    fprintf( stderr," Can not find any free markers. \"First_Free\" \n");
		    fprintf( stderr," Aborting\n");
		    exit(-1);
		}
		i = start;  /* restart going in opposite direction */
		step  = -step;
		i = i + step;
	    }
	}

	/* printf( " returns %i\n", i); */
	return i;
}

void process_input( )
{
	int i,c;
	double spaces[MAXLOC];

	scanf(" %i %i %i", &method.EvenCnt, &method.MidFlag, &method.CloseCnt);

/*	printf( " method:%i:%i:%i \n",method.EvenCnt, method.MidFlag, method.CloseCnt); */
	method.TotalCnt = method.EvenCnt*2 + method.MidFlag + method.CloseCnt;
   
	scanf(" %lf %lf %lf", &inc, &StartOffset, &StopOffset );
/*	printf( " inc=%f, Offsets=%f,%f\n", inc, StartOffset, StopOffset ); */

	scanf(" %1c ", &mapping );
	if ( mapping == 'k' ) mapping = 'K';
	if ( mapping == 'h' ) mapping = 'H';

	if ( (mapping != 'K') && (mapping != 'H') ) {
		fprintf( stderr," Mapping function not (K)osambi or (H)aldane\n");
		exit(-2);
	}

	i = 0;
	c = 1;

	while ( c == 1 ) {
	    c = scanf(" %lf", &spaces[i] );
	    if ( spaces[i] <= 0 ) {
		c = 0;
	    }
	    i++;
	    if ( i > MAXLOC ) {
		printf(" More markers than this program allows %i\n", MAXLOC );
		printf( " Aborting\n");
		exit( -1 );
	    } 

	}
	MarkerCnt = i;
	locations[0] = StartOffset;
	for ( i=0 ; i<MarkerCnt-1 ; i++ ) {
	   locations[i+1] = locations[i] + spaces[i];
	}
	locations[MarkerCnt] = locations[MarkerCnt-1] + StopOffset;

	if ( method.TotalCnt > MarkerCnt ) {
	    fprintf( stderr," Not Enough Markers to Impliment this Method.\n");
	    fprintf( stderr," Aborting\n");
	    exit( -1 );
	}

/*	printf(" marker locations:");
	for ( i=0 ; i<MarkerCnt ; i++ ) 
	    printf( " %6.3f",locations[i]);
	printf("\n"); */
}

void write_out( double loc ) {


	static int in_used[MAXLOC];
	int 	i, newmap, spot, j, cnt, flag;
	double	selloc[MAXLOC], foo;

	newmap = 0;  /* Scan to see if Markers used has changed, and off map run is needed*/
	for (i=0 ; i <= MarkerCnt ; i++ )
	    if (in_use[i] != in_used[i] ) newmap = 1;

	if ( newmap ) {  /* new to do a new/off map run */
	    printf("O       |" ); /* write marker map */
	    for ( i=0 ; i<=method.TotalCnt ;  i++ ) printf(" %3i", i+1 );
	    printf( "|"); 
                 /* Theta list */
	    printf( "  0.50000  ");
	    for ( i=0 ; i<MarkerCnt-1 ;  i++ ) {
	    	if ( (in_use[i] == USED) && (in_use[i+1] == USED) ) { 
		    printf( " %10.8f",to_theta(locations[i+1]-locations[i]) );
	        }
	    }
            printf( "| 0.50000  |");	 /* the foo value can be ignored */
		/* Now the markers to use */
	    for ( i=0 ; i<MarkerCnt ;  i++ ) {
	    	if ( in_use[i] == USED ) { 
		    printf( " %i ",i);
	        }
		in_used[i] = in_use[i];
	    }
	    printf( "\n");
	}


 /* write the regular on-map stuff */
	printf("%8.3f", loc );

/* Now to build a list for just in_use points */
	cnt = 0; /* build the list of just in_use loci */
	for ( i=0 ; i<MarkerCnt ; i++ ) {
           if ( in_use[i] == USED ) {
		selloc[cnt] = locations[i];
	        cnt++;
	   }
        }
/* write marker map */
	spot = -1;
	printf( "|" );
	for ( i=0 ; (i<cnt) ; i++ ) {
	   if ( (loc < selloc[i]) && (spot == -1) ) {
		spot = i;  /* spot now holds index loc is just infront of */
		printf( "   1");
	   }
	   printf( " %3i",i+2);
        }
	if ( spot == -1 ) {
	   printf( "   1");
	   spot = cnt;
        }
	printf( "|"); 

/* now add the trial loci to the array */
	if ( loc >= selloc[cnt-1] ) {
		selloc[cnt] = loc;
	        spot = cnt;
	} else {
	   spot = -1;
	   for ( i=0 ; (i<cnt)&&(spot==-1) ; i++ ) { 
	      if ( loc < selloc[i] ) {
	         for ( j=cnt ; j > i ; j-- ) { 
		    selloc[j] = selloc[j-1];
	         }
	         selloc[i] = loc;
	         spot = i;
	      }
	   }
	}
	for ( i=0 ; i<cnt ; i++ ) { /* and print out the selloc */
	   printf( " %10.8f", to_theta(selloc[i+1]-selloc[i]) );
	}
	printf( "|");
/* Now get and write the "Foo" value */
	if ( spot == 0 ) {
	   foo = selloc[spot+1] - selloc[spot];
	} else {
	   foo = selloc[spot] - selloc[spot-1];
	   if ( (foo == 0.0) && (spot < cnt) ) 
	      foo = (selloc[spot+1] - selloc[spot])/2.0;
	}
	printf( "%10.8f|", to_theta(foo) );	

/* Now write the selected Marker index values */
	for ( i=0 ; i<MarkerCnt ;  i++ ) {
	    if ( in_use[i] == USED ) { 
		printf( "%3i ",i);
	    }
	}
/*	printf( " :%i prev:%6.2f spot:%6.2f next:%6.2f", 
		spot, selloc[spot-1],selloc[spot],selloc[spot+1]); */
	printf( "\n");



}

void distribute( )
{
	int	i;
	double	loc, midpt, stophere;
	int	in_used[MAXLOC];
	int	flag;


	stophere = locations[MarkerCnt];

	for ( i =0 ; i < MAXLOC ; i++ )
		in_used[i] = -1;  /* should force offmap entry */

          /* The Big Loop, stepping through the locations */
	for ( loc = 0 ; loc <= stophere ; loc=loc+inc ) {

	   /* printf(" doing loc %f\n",loc); */
	    for ( i=0 ; i<MarkerCnt ; i++ ) 
		in_use[i] = FREE;


/*------------------ Do the Even Marker picking-------------------- */
	    if ( method.EvenCnt > 0 ) { 
	    	for ( i=1 ; i<=method.EvenCnt ; i++ ) {
		    /* Grab one on the left */
		    in_use[ First_Free( GO_LEFT, loc ) ] = 1;
		    /* Now grab one on the Right */
		    in_use[ First_Free( GO_RIGHT, loc ) ] = 1;
		}

	    }

/*------------------ Do the Midpoint Adding -----------------------*/
	    if ( method.MidFlag ) {
	    	i = where_at( loc );
	    	if ( i == 0 ) {
		    in_use[ First_Free( GO_RIGHT, loc ) ] = 1;
	    	} else if ( i == MarkerCnt ) {
		    in_use[ First_Free( GO_LEFT, loc ) ] = 1;
	    	} else {
		    midpt = locations[i-1]+((locations[i]-locations[i-1])/2.0); 
		    if ( loc > midpt ) {
		    	in_use[ First_Free( GO_RIGHT, loc ) ] = 1;
		    } else {
		    	in_use[ First_Free( GO_LEFT, loc ) ] = 1;
		    }
	    	}
	    }

/*----------Write out ----*/


	    write_out( loc );

	}

	printf( "EOF\n" );
}

int main( int argc, char *argv[] )
{
	process_input( );
	distribute();
}
        
