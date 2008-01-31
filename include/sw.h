#include <sys/resource.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define MAXSWNAME 32

struct swStopwatch {
  char swName[MAXSWNAME+1];
  struct rusage swStartRU;
  time_t swStartWallTime;
  struct rusage swAccumRU;
  time_t swAccumWallTime;
  int swStartedCount;
  int swRunning;
};

struct swStopwatch *swCreate(char *);
void swStart(struct swStopwatch *);
void swStop(struct swStopwatch *);
void swDump(struct swStopwatch *);
void swReset(struct swStopwatch *);
void *swMalloc(size_t size, char *fileName, int lineNo);
void swFree(void *pBlock, char *fileName, int lineNo);
void swDumpBlockUse();
