void dk_write2ptBRHeader (int loc1, int loc2);
void dk_write2ptBRData (double integral);
void dk_writeMPBRHeader ();
void dk_writeMPBRData (int posIdx, float traitPos, double ppl, double br);
void dk_writeMPMODHeader ();
void dk_writeMPMODData (int posIdx, float traitPos, double value, st_DKMaxModel *model);
void dk_write2ptMODHeader ();
void dk_write2ptMODData (char *description, double value, st_DKMaxModel *model);
void dk_copyMaxModel (double *arr, st_DKMaxModel *max);
