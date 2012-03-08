void initializeDB ();
void prepareDBStatements ();
long GetPedPosId (char *, int, double);
double GetDLikelihood (int, double, 
		double, double, double, double, 
		double, double, double, double, 
		double, double, double, double,
		int, int, double, int);
void SignOn (int, char *, int, char *);
void GetAnalysisId ();
void SignOff (int);
void SetDummyNullLikelihood ();
int CountWork (double, double);
int GetDWork (double, double, int, double *, char *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *);
void PutWork (int, double, int);

double GetQLikelihood (int, double, 
		       double, double, double, double,
		       double, double, double, double,
		       double, double, double, double,
		       double, double, double, double,
		       double, double, double, double,
		       double, double, double, double,
		       double, double, double,
		       int, int, double, int);
int GetQWork (double, double, int, double *, char *, double *, 
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *, double *,
	      double *, double *, double *);


double GetMarkerSetLikelihood(int pedPosId, int regionNo, int parentRegionNo, double parentRegionError, int parentRegionSplitDir);
int GetMarkerSetLikelihood_MCMC(int pedPosId);

