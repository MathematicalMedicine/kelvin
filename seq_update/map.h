#ifndef _MAP_H

#define BUFFLEN          256
#define MAX_MAP_CHR_LEN    8
#define MAX_MAP_NAME_LEN  32

#define MAP_CHR_COL        1
#define MAP_NAME_COL       2
#define MAP_ALTNAME_COL    3
#define MAP_AVGPOS_COL     4
#define MAP_MALEPOS_COL    5
#define MAP_FEMALEPOS_COL  6
#define MAP_BASEPAIR_COL   7


typedef struct {
  char chr[MAX_MAP_CHR_LEN],
    name[MAX_MAP_NAME_LEN],
    altname[MAX_MAP_NAME_LEN];
  double avgpos,
    malepos,
    femalepos;
  long basepair;
} st_mapmarker;

typedef struct {
  char chr[MAX_MAP_CHR_LEN];
  int nummarkers;
  st_mapmarker *markers;
} st_mapchr;

void read_map (char *filename);
st_mapmarker *find_mapmarker (char *chr, char *name);
st_mapmarker *next_mapmarker (st_mapmarker *mrk);
void dump_map ();
void free_map ();

#define _MAP_H
#endif
