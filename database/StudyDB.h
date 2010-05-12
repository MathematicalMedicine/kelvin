#include "mysql.h"

#define MAXSTMTLEN 256
#define MAXPARAMLEN 64

struct StudyDB {
  int studyId;
  int modelId;
  MYSQL *connection;
  char hostname[MAXPARAMLEN];
  char dBName[MAXPARAMLEN];
  char username[MAXPARAMLEN];
  char password[MAXPARAMLEN];
  char strAdhocStatement[MAXSTMTLEN];
  char strSelectModels[MAXSTMTLEN];
  char strInsertModels[MAXSTMTLEN];
  MYSQL_STMT *stmtSelectModels;
  MYSQL_STMT *stmtInsertModels;
  MYSQL_BIND bindModelsParams[5];
  MYSQL_BIND bindModelsResult[1];
  unsigned long dummyLength;
  my_bool dummyIsNull;
  my_bool dummyError;
  MYSQL_RES *resultSet;
  MYSQL_ROW row;
};

