/* Copyright (C) 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
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

#define MAX_DB_RETRIES 5

