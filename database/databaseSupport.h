void initializeDB ();
void prepareDBStatements ();
long GetPedPosId (char *, int, double);
double GetDLikelihood (int, double, 
		double, double, double, double, 
		double, double, double, double, 
		double, double, double, double,
		int, int, double, int);
void SignOn (int, char *, int, char *);
void SetDummyNullLikelihood ();
int CountWork (double, double);
int GetDWork (double, double, int, double *, char *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *);
void PutWork (int, double, int);

