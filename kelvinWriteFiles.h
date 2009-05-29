void writePPLFileHeader ();
void writePPLFileDetail (int dprime0Idx);
void write2ptBRFile (int loc1, int loc2);
void writeMPBRFileHeader ();
void writeMPBRFileDetail (int posIdx, float traitPos, float ppl, double avgLR);
void writeMPMODFileHeader ();
void writeMPMODFileDetail (int posIdx, float traitPos);
void writeMaximizingModel (char *modelDescription, double myMOD, int myDPrimeIdx,
			   int myThetaIdx);
void write2ptMODFile (int loc1, int loc2, int dprime0Idx);







