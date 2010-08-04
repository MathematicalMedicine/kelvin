#include "mysql.h"

#define MAXSTMTLEN 512
#define MAXPARAMLEN 64

struct StudyDB {
  // Common...
  char role[MAXPARAMLEN];
  int studyId;
  char hostname[MAXPARAMLEN];
  char dBName[MAXPARAMLEN];
  char username[MAXPARAMLEN];
  char password[MAXPARAMLEN];
  char pedigreeRegEx[33];
  MYSQL *connection;
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
  // GetDLOD...
  MYSQL_STMT *stmtGetDLOD;
  MYSQL_BIND bindGetDLOD[15];
  char strGetDLOD[MAXSTMTLEN];
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
  char algorithm[3];
  char programVersion[33];
  int serverId;
  // GetWork...
  MYSQL_STMT *stmtGetWork;
  MYSQL_BIND bindGetWork[3];
  char strGetWork[MAXSTMTLEN];
  double lowPosition;
  double highPosition;
  // GetDParts...
  MYSQL_STMT *stmtGetDParts;
  char strGetDParts[MAXSTMTLEN];
  // GetWorkResults
  MYSQL_STMT *stmtGetWorkResults;
  MYSQL_BIND bindGetWorkResults[16];
  char strGetWorkResults[MAXSTMTLEN];
  double pedTraitPosCM;
  // PutWork...
  MYSQL_STMT *stmtPutWork;
  MYSQL_BIND bindPutWork[3];
  char strPutWork[MAXSTMTLEN];
  // other...
  int bogusLODs;
  int realLODs;
};

