#include "mysql.h"

struct StudyDB {
  int studyId;
  char hostname[64];
  char dBName[64];
  char username[64];
  char password[64];
  char statement[256];
  MYSQL *connection;
  MYSQL_STMT *parsedStatement;
  MYSQL_RES *resultSet;
  MYSQL_ROW row;
};

