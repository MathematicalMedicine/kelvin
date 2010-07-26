#include "mysql.h"

#define MAXSTMTLEN 256
#define MAXPARAMLEN 64

struct StudyDB {
  // Common...
  char role[MAXPARAMLEN];
  int studyId;
  char hostname[MAXPARAMLEN];
  char dBName[MAXPARAMLEN];
  char username[MAXPARAMLEN];
  char password[MAXPARAMLEN];
  MYSQL *connection;
  // Adhoc...
  char strAdhocStatement[MAXSTMTLEN];
  MYSQL_RES *resultSet;
  MYSQL_ROW row;
  // GetPedPosId...
  MYSQL_STMT *stmtGetPedPosId;
  MYSQL_BIND bindGetPedPosId[4];
  char strGetPedPosId[MAXSTMTLEN];
  char pedigreeSId[16];
  int chromosomeNo;
  double inRefTraitPosCM;
  MYSQL_BIND bindGetPedPosIdResults[1];
  int pedPosId;
  // GetDLOD...
  MYSQL_STMT *stmtGetDLOD;
  MYSQL_BIND bindGetDLOD[15];
  char strGetDLOD[MAXSTMTLEN];
  double inDGF;
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
  // GetDLOD results...
  MYSQL_STMT *stmtGetDLODResults;
  MYSQL_BIND bindGetDLODResults[3];
  char strGetDLODResults[MAXSTMTLEN];
  int regionId;
  int markerCount;
  double lOD;
  // Sign-on...
  MYSQL_STMT *stmtSignOn;
  MYSQL_BIND bindSignOn[6];
  char strSignOn[MAXSTMTLEN];
  char pedigreeRegEx[32];
  char algorithm[2];
  char programVersion[32];
  // GetWork...
  MYSQL_STMT *stmtGetWork;
  MYSQL_BIND bindGetWork[3];
  char strGetWork[MAXSTMTLEN];
  int serverId;
  double lowPosition;
  double highPosition;
  // GetDTParts...
  MYSQL_STMT *stmtGetDTParts;
  char strGetDTParts[MAXSTMTLEN];
  // GetWorkResults
  MYSQL_STMT *stmtGetWorkResults;
  MYSQL_BIND bindGetWorkResults[16];
  char strGetWorkResults[MAXSTMTLEN];
  double outPedTraitPosCM;
  double outDGF;
  // PutWork...
  MYSQL_STMT *stmtPutWork;
  MYSQL_BIND bindPutWork[2];
  char strPutWork[MAXSTMTLEN];
  // other...
  int bogusLODs;
  int realLODs;
};

