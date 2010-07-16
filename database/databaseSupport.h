void initializeDB ();
void prepareDBStatements ();
long GetPedPosId (char *inPedigreeSId, int inChromosomeNo, double inRefTraitPosCM);
double GetDLOD (int inPedPosId, double inDGF,
		double inLC1BigPen, double inLC1BigLittlePen, double inLC1LittleBigPen, double inLC1LittlePen,
		double inLC2BigPen, double inLC2BigLittlePen, double inLC2LittleBigPen, double inLC2LittlePen,
		double inLC3BigPen, double inLC3BigLittlePen, double inLC3LittleBigPen, double inLC3LittlePen,
		int inRegionId);
