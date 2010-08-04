void initializeDB ();
void prepareDBStatements ();
long GetPedPosId (char *, int, double);
double GetDLOD (int, double, 
		double, double, double, double, 
		double, double, double, double, 
		double, double, double, double,
		int);
void SignOn (int, char *, int, char *);
int GetDWork (double, double, double *, char *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *);
void PutWork (int, double);

