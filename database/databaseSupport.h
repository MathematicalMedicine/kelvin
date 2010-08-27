void initializeDB ();
void prepareDBStatements ();
long GetPedPosId (char *, int, double);
double GetDAltL (int, double, 
		double, double, double, double, 
		double, double, double, double, 
		double, double, double, double,
		int, int, double, int);
void SignOn (int, char *, int, char *);
int CountWork (double, double);
int GetDWork (double, double, double *, char *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *);
void PutWork (int, double, int);

