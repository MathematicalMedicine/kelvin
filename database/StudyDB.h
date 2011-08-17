#include "mysql.h"

#define MAXSTMTLEN 512
#define MAXPARAMLEN 64

struct StudyDB {
  // Common...
  char role[MAXPARAMLEN];
  int studyId;
  char dBHostname[MAXPARAMLEN];
  char dBName[MAXPARAMLEN];
  char username[MAXPARAMLEN];
  char password[MAXPARAMLEN];
  char pedigreeRegEx[33];
  char pedigreeNotRegEx[33];
  MYSQL *connection;
  int driverPosIdx;
  // Adhoc...
  char strAdhocStatement[MAXSTMTLEN];
  MYSQL_RES *resultSet;
  MYSQL_ROW row;
  // GetPedPosId...
  MYSQL_STMT *stmtGetPedPosId;
  MYSQL_BIND bindGetPedPosId[4];
  char strGetPedPosId[MAXSTMTLEN];
  char pedigreeSId[17];
  int chromosomeNo;
  double refTraitPosCM;
  MYSQL_BIND bindGetPedPosIdResults[1];
  int pedPosId;
  // GetDLikelihood...
  int posEvals;
  MYSQL_STMT *stmtGetDLikelihood;
  MYSQL_BIND bindGetDLikelihood[18];
  char strGetDLikelihood[MAXSTMTLEN];
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
  char hostName[33];
  int processId;
  int keepAliveFlag;
  char algorithm[3];
  char programVersion[33];
  int serverId;
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
  MYSQL_BIND bindCountWorkResults[16];
  char strCountWorkResults[MAXSTMTLEN];
  long workCount;
  // GetWork...
  MYSQL_STMT *stmtGetWork;
  MYSQL_BIND bindGetWork[4];
  char strGetWork[MAXSTMTLEN];
  double lowPosition;
  double highPosition;
  int locusListType;
  // GetDParts...
  MYSQL_STMT *stmtGetDParts;
  char strGetDParts[MAXSTMTLEN];
  // GetDWorkResults
  MYSQL_STMT *stmtGetDWorkResults;
  MYSQL_BIND bindGetDWorkResults[16];
  char strGetDWorkResults[MAXSTMTLEN];
  double pedTraitPosCM;
  // PutWork...
  MYSQL_STMT *stmtPutWork;
  MYSQL_BIND bindPutWork[4];
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
  MYSQL_BIND bindGetQLikelihood[33];
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

  // GetQWorkResults
  MYSQL_STMT *stmtGetQWorkResults;
  MYSQL_BIND bindGetQWorkResults[31];
  char strGetQWorkResults[MAXSTMTLEN];
  //double pedTraitPosCM;
};

