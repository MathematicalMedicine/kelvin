
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

void dump_lDT (double **);
int saveTrait (int, char *, double **);
int restoreTrait (int, char *, double **);
int saveMarker (char *, int, int, char **, double *);
int restoreMarker (char *, int, int, char **, double *);
int saveAlternative (char *, int, double, double **);
int restoreAlternative (char *, int, double, double **);
