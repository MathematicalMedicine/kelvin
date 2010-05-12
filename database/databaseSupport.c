#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "../utils/sw.h"
#include "../utils/utils.h"

#include "StudyDB.h"

typedef struct
{
  double DD, Dd, dD, dd, DDSD, DdSD, dDSD, ddSD, threshold;
} st_DKMaxModelPenVector;

typedef struct
{
  int posIdx;   // which stores loc2 for 2pt and posIdx for mp
  double *dprime, theta[2], alpha, dgf, mf, r2;
  st_DKMaxModelPenVector *pen;
} st_DKMaxModel;

extern st_DKMaxModel dk_curModel;
extern struct StudyDB studyDB;

void initializeDB () {

  /* Initialize structure for a MySQL connection. */
  if ((studyDB.connection = mysql_init(NULL)) == NULL)
    FATAL("Cannot initialize MySQL (%s)", strerror(errno));

  /* Connect. */
  if (!mysql_real_connect(studyDB.connection, studyDB.hostname, studyDB.username, studyDB.password, NULL, 0, NULL, 0))
    ERROR("Cannot connect to MySQL on hostname [%s] as username [%s/%s] (%s)", studyDB.hostname, 
	  studyDB.username, studyDB.password, mysql_error(studyDB.connection));

  /* Change database. */
  if (mysql_select_db(studyDB.connection, studyDB.dBName))
    ERROR("Cannot change MySQL db (%s)", mysql_error(studyDB.connection));

  /* Verify our studyId. */
  sprintf (studyDB.strAdhocStatement, "Select Description from Studies where StudyId = %d", studyDB.studyId);
  if (mysql_query (studyDB.connection, studyDB.strAdhocStatement))
    ERROR("Cannot select study information (%s:%s)", studyDB.strAdhocStatement, mysql_error(studyDB.connection));
  if ((studyDB.resultSet = mysql_store_result (studyDB.connection)) == NULL)
    ERROR("Cannot retrieve study information (%s)", mysql_error(studyDB.connection));
  if (mysql_num_rows (studyDB.resultSet) == 0)
    ERROR("Study %d not found", studyDB.studyId);
  else {
    if ((studyDB.row = mysql_fetch_row (studyDB.resultSet)) == NULL)
      ERROR("Cannot fetch study information (%s)", mysql_error(studyDB.connection));
    INFO ("Storing/retrieving results under study %d (%s)", studyDB.studyId, studyDB.row[0]);
  }
}

void *prepareModels () {

  // Prepare the Models select and insert statements
  studyDB.stmtInsertModels = mysql_stmt_init (studyDB.connection);
  studyDB.stmtSelectModels = mysql_stmt_init (studyDB.connection);

  // These two statements share parameter bindings
  memset (studyDB.bindModelsParams, 0, sizeof(studyDB.bindModelsParams));

  studyDB.bindModelsParams[0].buffer_type = MYSQL_TYPE_DOUBLE;
  studyDB.bindModelsParams[0].buffer = &dk_curModel.alpha;
  studyDB.bindModelsParams[0].buffer_length = sizeof (dk_curModel.alpha);
  studyDB.bindModelsParams[0].length = 0; // Ignored for numeric and temporal data
  studyDB.bindModelsParams[0].is_null = (my_bool *) 0; // Not null - notice this flag is a pointer!
  studyDB.bindModelsParams[0].is_unsigned = 0; // Not unsigned (client-side)

  studyDB.bindModelsParams[1].buffer_type = MYSQL_TYPE_DOUBLE;
  studyDB.bindModelsParams[1].buffer = &dk_curModel.dgf;
  studyDB.bindModelsParams[1].buffer_length = sizeof (dk_curModel.dgf);
  studyDB.bindModelsParams[1].length = 0; // Ignored for numeric and temporal data
  studyDB.bindModelsParams[1].is_null = (my_bool *) 0; // Not null - notice this flag is a pointer!
  studyDB.bindModelsParams[1].is_unsigned = 0; // Not unsigned (client-side)

  studyDB.bindModelsParams[2].buffer_type = MYSQL_TYPE_DOUBLE;
  studyDB.bindModelsParams[2].buffer = &dk_curModel.pen->DD;
  studyDB.bindModelsParams[2].buffer_length = sizeof (dk_curModel.pen->DD);
  studyDB.bindModelsParams[2].length = 0; // Ignored for numeric and temporal data
  studyDB.bindModelsParams[2].is_null = (my_bool *) 0; // Not null - notice this flag is a pointer!
  studyDB.bindModelsParams[2].is_unsigned = 0; // Not unsigned (client-side)

  studyDB.bindModelsParams[3].buffer_type = MYSQL_TYPE_DOUBLE;
  studyDB.bindModelsParams[3].buffer = &dk_curModel.pen->Dd;
  studyDB.bindModelsParams[3].buffer_length = sizeof (dk_curModel.pen->Dd);
  studyDB.bindModelsParams[3].length = 0; // Ignored for numeric and temporal data
  studyDB.bindModelsParams[3].is_null = (my_bool *) 0; // Not null - notice this flag is a pointer!
  studyDB.bindModelsParams[3].is_unsigned = 0; // Not unsigned (client-side)

  studyDB.bindModelsParams[4].buffer_type = MYSQL_TYPE_DOUBLE;
  studyDB.bindModelsParams[4].buffer = &dk_curModel.pen->dd;
  studyDB.bindModelsParams[4].buffer_length = sizeof (dk_curModel.pen->dd);
  studyDB.bindModelsParams[4].length = 0; // Ignored for numeric and temporal data
  studyDB.bindModelsParams[4].is_null = (my_bool *) 0; // Not null - notice this flag is a pointer!
  studyDB.bindModelsParams[4].is_unsigned = 0; // Not unsigned (client-side)

  strncpy (studyDB.strSelectModels, "Select ModelId from Models where Alpha = ? AND DGF = ? AND BigPen = ? AND BigLittlePen = ? AND LittlePen = ?", MAXSTMTLEN-1);
  strncpy (studyDB.strInsertModels, "Insert into Models (Alpha, DGF, BigPen, BigLittlePen, LittlePen) values (?,?,?,?,?)", MAXSTMTLEN-1);

  if (mysql_stmt_prepare (studyDB.stmtSelectModels, studyDB.strSelectModels, strlen (studyDB.strSelectModels)))
    ERROR("Cannot prepare Models selection statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_prepare (studyDB.stmtInsertModels, studyDB.strInsertModels, strlen (studyDB.strInsertModels)))
    ERROR("Cannot prepare Models insertion statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtSelectModels, studyDB.bindModelsParams))
    ERROR("Cannot bind Models selection statement (%s)", mysql_error(studyDB.connection));
  if (mysql_stmt_bind_param (studyDB.stmtInsertModels, studyDB.bindModelsParams))
    ERROR("Cannot bind Models insertion statement (%s)", mysql_error(studyDB.connection));

  memset (studyDB.bindModelsResult, 0, sizeof(studyDB.bindModelsResult));
  studyDB.bindModelsResult[0].buffer_type = MYSQL_TYPE_LONG;
  studyDB.bindModelsResult[0].buffer = &studyDB.modelId;

  if (mysql_stmt_bind_result (studyDB.stmtSelectModels, studyDB.bindModelsResult))
    ERROR("Cannot bind result for Models select statement (%s)", mysql_error(studyDB.connection));

}

void getCurrentModelId () {
  
  // Get the ModelId matching the current set of trait parameters
  if (mysql_stmt_execute (studyDB.stmtSelectModels))
    ERROR("Cannot execute Models select statement (%s)", mysql_stmt_error(studyDB.stmtSelectModels));
  if (mysql_stmt_store_result (studyDB.stmtSelectModels))
    ERROR("Cannot retrieve Model information (%s)", mysql_stmt_error(studyDB.stmtSelectModels));

  if (!mysql_stmt_fetch (studyDB.stmtSelectModels))
    fprintf (stderr, "%d ", studyDB.modelId);
  else {
    // Didn't find it, do an insert
    if (mysql_stmt_execute (studyDB.stmtInsertModels))
      ERROR("Cannot execute Models insertion statement (%s)", mysql_stmt_error(studyDB.stmtInsertModels));
    if (mysql_stmt_execute (studyDB.stmtSelectModels))
      ERROR("Cannot execute Models select statement (%s)", mysql_stmt_error(studyDB.stmtInsertModels));
    if (mysql_stmt_store_result (studyDB.stmtSelectModels))
      ERROR("Cannot retrieve Model information (%s)", mysql_stmt_error(studyDB.stmtInsertModels));
    if (!mysql_stmt_fetch (studyDB.stmtSelectModels))
      fprintf (stderr, "+%d ", studyDB.modelId);
    else
      ERROR("Cannot fetch Model information even after insert (%s)", mysql_stmt_error(studyDB.stmtSelectModels));
  }
}
