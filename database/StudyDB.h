#include "mysql.h"

#define MAXSTMTLEN 256
#define MAXPARAMLEN 64

struct StudyDB {
  // Common...
  char role[MAXPARAMLEN];
  int inStudyId;
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
  char inPedigreeSId[16];
  int inChromosomeNo;
  double inRefTraitPosCM;
  MYSQL_BIND bindGetPedPosIdResults[1];
  double outPedPosId;
  // GetDLOD...
  MYSQL_STMT *stmtGetDLOD;
  MYSQL_BIND bindGetDLOD[16];
  char strGetDLOD[MAXSTMTLEN];
  int inPedPosId;
  double inAlpha;
  double inDGF;
  double inLC1BigPen;
  double inLC1BigLittlePen;
  double inLC1LittleBigPen;
  double inLC1LittlePen;
  double inLC2BigPen;
  double inLC2BigLittlePen;
  double inLC2LittleBigPen;
  double inLC2LittlePen;
  double inLC3BigPen;
  double inLC3BigLittlePen;
  double inLC3LittleBigPen;
  double inLC3LittlePen;
  int inRegionId;
  MYSQL_STMT *stmtGetDLODResults;
  MYSQL_BIND bindGetDLODResults[2];
  char strGetDLODResults[MAXSTMTLEN];
  int outRegionId;
  double outLOD;
};

