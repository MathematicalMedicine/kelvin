
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#ifndef __TOOLS_H__
#define __TOOLS_H__

/* number of bits in an integer */
#define INT_BITS		(sizeof(int)*8)

int is_line_blank_or_comment (char *line);
int is_line_blank (char *line);
int is_line_comment (char *line);
char *get_nonblank_line (char *pLine, int maxLen, FILE * fp, int *pLineNo);
char *fgetlongs (char **buff, int *bufflen, FILE * fp);


#endif
