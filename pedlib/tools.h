#ifndef __TOOLS_H__
#define __TOOLS_H__

/* number of bits in an integer */
#define INT_BITS		(sizeof(int)*8)

int is_line_blank_or_comment (char *line);
int is_line_blank (char *line);
int is_line_comment (char *line);
char *get_nonblank_line (char *pLine, int maxLen, FILE * fp, int *pLineNo);

#endif
