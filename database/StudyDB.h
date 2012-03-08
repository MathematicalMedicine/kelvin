#include "mysql.h"
#include <regex.h>

#define MAXSTMTLEN 512
#define MAXPARAMLEN 64

struct StudyDB {
  // Common...
  char role[MAXPARAMLEN];
  char studyLabel[MAXPARAMLEN];
  char studyDescription[128];
  int studyId;
  int liabilityClassCnt;
  char imprintingFlag[2];
  char dBHostname[MAXPARAMLEN];
  char dBName[MAXPARAMLEN];
  char username[MAXPARAMLEN];
  char password[MAXPARAMLEN];
  char pedigreeRegEx[33];
  char pedigreeNotRegEx[33];
  int analysisId;
  regex_t includePattern;
  regex_t excludePattern;
  regmatch_t pmatch[1];
  MYSQL *connection;
  int driverPosIdx;
  // Adhoc...
  char strAdhocStatement[MAXSTMTLEN];
  MYSQL_RES *resultSet;
  MYSQL_ROW row;
  // GetStudyId...
  MYSQL_STMT *stmtGetStudyId;
  MYSQL_BIND bindGetStudyId[4];
  char strGetStudyId[MAXSTMTLEN];
  // GetAnalysisId...
  MYSQL_STMT *stmtGetAnalysisId;
  MYSQL_BIND bindGetAnalysisId[4];
  char strGetAnalysisId[MAXSTMTLEN];

  // GetPedPosId...
  MYSQL_STMT *stmtGetPedPosId;
  MYSQL_BIND bindGetPedPosId[5];
  char strGetPedPosId[MAXSTMTLEN];
  char pedigreeSId[17];
  int chromosomeNo;
  double refTraitPosCM;
  double refTraitPosCMleft;
  double refTraitPosCMright;
  MYSQL_BIND bindGetPedPosIdResults[1];
  int pedPosId;
  int posEvals;
  // allocate space to store markerset likelihood for each MCMC sampling
  // unfortunately we have to keep them to calculate the LR per sample, 
  // otherwise we get wrong results
  // just server side though
  int markerSetLikelihoodFlag;
  int markerSetPedPosId;
  double *markerSetLikelihood;
  char strGetMarkerSetLikelihood_MCMC[MAXSTMTLEN];
  MYSQL_STMT *stmtGetMarkerSetLikelihood_MCMC;
  char strGetMarkerSetId[MAXSTMTLEN];
  MYSQL_STMT *stmtGetMarkerSetId;
  MYSQL_BIND bindGetMarkerSetId[1];
  // GetMarkerSetLiklihood
  MYSQL_STMT *stmtGetMarkerSetLikelihood;
  MYSQL_BIND bindGetMarkerSetLikelihood[6];
  char strGetMarkerSetLikelihood[MAXSTMTLEN];
  // GetMarkerSetLikelihood results
  MYSQL_STMT *stmtGetMarkerSetLikelihoodResults;
  MYSQL_BIND bindGetMarkerSetLikelihoodResults[3];
  char strGetMarkerSetLikelihoodResults[MAXSTMTLEN];
  
  // GetDLikelihood...
  MYSQL_STMT *stmtGetDLikelihood;
  MYSQL_BIND bindGetDLikelihood[19];
  char strGetDLikelihood[MAXSTMTLEN];
  // index - LC
  int partsIdx[3];
  double dGF;
  double lC1BigPen;
  double lC1BigLittlePen;
  double lC1LittleBigPen;
  double lC1LittlePen;
  double lC2BigPen;
  double lC2BigLittlePen;
  double lC2LittleBigPen;
  double lC2LittlePen;
  double lC3BigPen;
  double lC3BigLittlePen;
  double lC3LittleBigPen;
  double lC3LittlePen;
  int regionNo;
  int parentRegionNo;
  double parentRegionError;
  int parentRegionSplitDir;
  // GetDLikelihood results...
  MYSQL_STMT *stmtGetDLikelihoodResults;
  MYSQL_BIND bindGetDLikelihoodResults[3];
  char strGetDLikelihoodResults[MAXSTMTLEN];
  int regionId;
  int markerCount;
  double lOD;
  // Sign-on...
  MYSQL_STMT *stmtSignOn;
  MYSQL_BIND bindSignOn[12];
  char strSignOn[MAXSTMTLEN];
  MYSQL_BIND bindSignOnResults[1];
  MYSQL_STMT *stmtSignOnResults;
  char strSignOnResults[MAXSTMTLEN];
  char hostName[33];
  int processId;
  int keepAliveFlag;
  char algorithm[3];
  char programVersion[33];
  int serverId;
  // Sign-off...
  MYSQL_STMT *stmtSignOff;
  MYSQL_BIND bindSignOff[2];
  char strSignOff[MAXSTMTLEN];
  int exitStatus;
  // SetDummyNullLikelihood...
  MYSQL_STMT *stmtSetDummyNullLikelihood;
  MYSQL_BIND bindSetDummyNullLikelihood[1];
  char strSetDummyNullLikelihood[MAXSTMTLEN];
  // CountWork...
  MYSQL_STMT *stmtCountWork;
  MYSQL_BIND bindCountWork[3];
  char strCountWork[MAXSTMTLEN];
  // CountWorkResults
  MYSQL_STMT *stmtCountWorkResults;
  MYSQL_BIND bindCountWorkResults[1];
  char strCountWorkResults[MAXSTMTLEN];
  long workCount;
  // GetWork...
  MYSQL_STMT *stmtGetWork;
  MYSQL_BIND bindGetWork[7];
  char strGetWork[MAXSTMTLEN];
  double lowPosition;
  double highPosition;
  int locusListType;
  MYSQL_BIND bindGetWorkResults[6];
  // GetDParts...
  MYSQL_STMT *stmtGetDParts;
  MYSQL_BIND bindGetDParts[13];
  MYSQL_BIND bindGetDPartsResults[13];
  char strGetDParts[MAXSTMTLEN];
  // GetDWorkResults
  MYSQL_STMT *stmtGetDWorkResults;
  MYSQL_BIND bindGetDWorkResults[16];
  char strGetDWorkResults[MAXSTMTLEN];
  double pedTraitPosCM;
  // PutWork...
  MYSQL_STMT *stmtPutWork;
  MYSQL_BIND bindPutWork[5];
  char strPutWork[MAXSTMTLEN];
  int runtimeCostSec;
  // other...
  int bogusLikelihoods;
  int realLikelihoods;

  // MCMC
  int totalSampleCount;
  int sampleIdStart;
  int sampleIdEnd;
  int MCMC_flag;

  int traitType;

  // GetQLikelihood...
  //  int posEvals;
  MYSQL_STMT *stmtGetQLikelihood;
  MYSQL_BIND bindGetQLikelihood[34];
  char strGetQLikelihood[MAXSTMTLEN];
  //  double dGF;
  double lC1BigMean;
  double lC1BigLittleMean;
  double lC1LittleBigMean;
  double lC1LittleMean;
  double lC2BigMean;
  double lC2BigLittleMean;
  double lC2LittleBigMean;
  double lC2LittleMean;
  double lC3BigMean;
  double lC3BigLittleMean;
  double lC3LittleBigMean;
  double lC3LittleMean;
  double lC1BigSD;
  double lC1BigLittleSD;
  double lC1LittleBigSD;
  double lC1LittleSD;
  double lC2BigSD;
  double lC2BigLittleSD;
  double lC2LittleBigSD;
  double lC2LittleSD;
  double lC3BigSD;
  double lC3BigLittleSD;
  double lC3LittleBigSD;
  double lC3LittleSD;
  double lC1Threshold;
  double lC2Threshold;
  double lC3Threshold;
  //int regionNo;
  //int parentRegionNo;
  //double parentRegionError;
  //int parentRegionSplitDir;

  // GetQLikelihood results...
  MYSQL_STMT *stmtGetQLikelihoodResults;
  MYSQL_BIND bindGetQLikelihoodResults[3];
  char strGetQLikelihoodResults[MAXSTMTLEN];
  //int regionId;
  //int markerCount;
  //double lOD;

  // GetQParts...
  MYSQL_STMT *stmtGetQParts;
  char strGetQParts[MAXSTMTLEN];
  MYSQL_BIND bindGetQParts[28];
  MYSQL_BIND bindGetQPartsResults[28];

  // GetQWorkResults
  MYSQL_STMT *stmtGetQWorkResults;
  MYSQL_BIND bindGetQWorkResults[31];
  char strGetQWorkResults[MAXSTMTLEN];
  //double pedTraitPosCM;
};

