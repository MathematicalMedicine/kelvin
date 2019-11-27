/**********************************************************************
 * NICE Utilities
 * Alberto Maria Segre
 * General Definitions, Software Timers and Message Logging System.
 * 
 * Copyright (c) 1996-2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "niceaux.h"
#include <time.h>		/* Needed for seeding and timing */
#include <stdlib.h>		/* Needed for malloc and random */
#include <signal.h>		/* Needed for signal */
#include <unistd.h>		/* Needed for alarm */
#include <stdio.h>		/* Needed for fprintf */
#include <stdarg.h>
#ifndef WIN32
#include <syslog.h>
#endif

#ifndef NOALARM
/**********************************************************************
 * Software timers provide a way for the NICE infrastructure (or any
 * application using this library) to multiplex multiple signal
 * handlers on a single signal (say this three times fast). This is
 * useful because each Unix process has only one alarm available.
 *
 * Two types of software timers are supported by the alarm mechanism:
 * an interrupt, which invokes a handler function and then continues
 * execution (after optionally rescheudling the interrupt), and a
 * timeout, which aborts the current thread of execution and unwinds
 * back to where the timeout was set (nonlocal exit).
 *
 * Timers and their handlers are maintained in temporal order on a
 * linked list of cells. The lone Unix alarm signal is used to trigger
 * processing of the first handler on the linked list, and then the
 * alarm is reset as appropriate.
 **********************************************************************/
#define MINALARMINT	1	/* Granularity of alarm system, secs */
#define NONLOCAL	1	/* 00000001 */
#define RECURRING	2	/* 00000010 */

/**********************************************************************
 * We'll use a linked list to store pointers to pending alarms in
 * temporal order. A better solution would be to use a heap or
 * priority queue, or even, for example, Berkeley DB's internal BTree
 * data structure; but since there are probably going to be only a few
 * alarms set at any one time O(N) is good enough.
 **********************************************************************/
typedef struct alarmCell
{
  int id;		/* Alarm ID. */
  alarmhandler handler;	/* Alarm handler. */
  unsigned int props;	/* Alarm properties. */
  unsigned int time;	/* Alarm timeout. */

  unsigned int delay;	/* Time to next alarm. */
  struct alarmCell *next;
}
AlarmCell;

/* Alarm variables. */
AlarmCell *alarms = NULL;	/* List of alarm pointers. */
AlarmCell *freeAlarms = NULL;	/* Free list of alarm pointers. */
unsigned int alarmIndex = 0;	/* Alarm id counter. */

/**********************************************************************
 * There is one subtle problem here. Since we suspend the alarms while
 * manipulating the alarm list, and since the granularity of the timer
 * is on the order of a second, if you add/remove a lot of alarms, it
 * is possible that the timer might in effect never advance, being
 * unset after a small fraction of a second, and restarted at the
 * whole section shortly thereafter. Everything works fine if you're
 * handling only a few alarms a minute, but applications with heavy
 * communication needs (e.g., if each socket primitive were to spawn a
 * new alarm, leading to many alarms per second) may suffer from this
 * problem.
 *
 * Our solution is to maintain an estimate of the rate of alarm
 * manipulations per second, and then correct, probabilistically, for
 * this factor when setting/clearing an alarm. We'll get an estimate
 * by dividing number of alarm operations by elapsed time.
 **********************************************************************/
unsigned int alarmOpRate = 20;	/* Number of alarm ops per second. */
unsigned int alarmOpInt = 100;  /* Alarm estimation interval. */
unsigned int alarmOpCount = 0;	/* Current count of alarm ops. */
time_t alarmOpLastTime;		/* Used to estimate alarm op rate. */

/* Alarm function prototypes for internal use only. */
int niceSetAlarm (unsigned int time, alarmhandler handler, 
		  int recurs, int returns);
unsigned int niceAddCell (AlarmCell *new, unsigned int timeLeft);
int niceRemoveCell (int id, int timeLeft);
static void niceFireAlarm (int signum);
void niceAlarmEstimate ();

/* Environment and default handler for nonlocal exit alarms (e.g.,
 * timeouts). */
static jmp_buf timeoutenv;
static void niceNonLocalAlarmHandler (int);

/**********************************************************************
 * Add a new alarm, returning an alarm identifier (a small
 * integer). The identifier can be used to cancel the pending alarm.
 *
 * The alarm handler almost follows sigaction() semantics; the only
 * difference between alarmhandler and sighandler_t is that the signal
 * number is unncessary -- it would always be the same.
 *
 * If recurs is TRUE, the alarm, upon firing, is immediately
 * rescheduled for the original time interval. Else the alarm is a
 * one-shot deal.
 *
 * This function is for internal use only. It is called by both
 * niceSetInterrupt () and niceSetTimeout ().
 **********************************************************************/
int
niceSetAlarm (unsigned int time, alarmhandler handler, int recurs, int returns)
{
  AlarmCell *new;
  struct sigaction sigact;			/* Signals. */
#ifdef TMRDBG
  unsigned int tmp;
  AlarmCell *cell;
#endif
  /* Turn off alarm handling while manipulating the alarms
   * list. Remember how much time was left to the first alarm so you
   * can restart things once you've handled the insertion. */
  unsigned int timeLeft = alarm(0);

  /* Sanity check: no pending alarms and timeLeft is not zero. Someone
   * else is fooling around with your timer! */
  if (timeLeft != 0 && !alarms)
    {
#ifdef TMRDBG
      fprintf (stderr, "Error: corrupted timer.\n");
#endif
      return (ERROR);
    }

  /* If there are no other alarms already pending, set up the universal
   * signal handler. */
  if (!alarms)
    {
      sigact.sa_handler = niceFireAlarm;
      sigemptyset (&sigact.sa_mask);
      sigact.sa_flags = 0;
      if (sigaction (SIGALRM, &sigact, NULL) == ERROR)
	{
#ifdef TMRDBG
	  fprintf (stderr, "Error: can't set signal handler.\n");
#endif
	  return (ERROR);
	}
    }

  /* Update your rate estimate, if necessary. */
  if (++alarmOpCount > alarmOpInt)
    niceAlarmEstimate ();

  /* Create (or recyle) a new cell for the alarm. */
  if (!freeAlarms)
    {
      new = (AlarmCell *) malloc (sizeof(AlarmCell));
      new->next = NULL;
    }
  else
    {
      new = freeAlarms;
      freeAlarms = new->next;
      new->next = NULL;
    }
  /* OK, set up the appropriate contents in the new cell. */
  new->id = ++alarmIndex;
  new->handler = handler;
  new->props = ((0 | (returns?0:NONLOCAL)) | (recurs?RECURRING:0));
  new->time = MAX(time, MINALARMINT);

  /* Insert the cell in the appropriate place, and update time to
   * first firing. */
  timeLeft = niceAddCell (new, timeLeft);

  /* Make estimated alarm rate correction. */
  if (timeLeft > 0 && !RANDINT(alarmOpRate))
    {
#if FALSE
      fprintf (stderr, "Alarm: adjusting.\n");
#endif      
      timeLeft--;
    }
    
#ifdef TMRDBG
  tmp = timeLeft;
  fprintf (stderr, "Alarm [+%d]:", new->id);
  cell = alarms;
  while (cell)
    {
      fprintf (stderr, " [%d@%u:%u]%s%s", cell->id, tmp, cell->delay,
	       (cell->props & RECURRING)?"*":"", (cell->props & NONLOCAL)?"!":"");
      tmp = tmp + cell->delay;
      cell = cell->next;
    }
  fprintf (stderr, "\n");
#endif

  /* If there are alarms ready to fire, do so (TODO: what happens if
   * the alarm that fires is nonreturning? Is that a problem? Or
   * not?). Otherwise, don't forget to reset the alarm before you
   * exit. */
  if (timeLeft == 0)
    niceFireAlarm (0);
  else
    alarm (timeLeft);
  return (new->id);
}

/**********************************************************************
 * Set an interrupt. Internally, this just calls niceSetAlarm () to
 * establish a normally-returning alarm handler. Returns a positive
 * integer alarm number, which can be used to cancel the pending
 * interrupt.
 **********************************************************************/
int
niceSetInterrupt (unsigned int time, alarmhandler handler, int recurs)
{
  return (niceSetAlarm (time, handler, recurs, TRUE));
}

/**********************************************************************
 * Set a timeout. Internally, this just calls niceSetAlarm () to
 * establish a non-locally returning alarm handler using the
 * appropriate alarm handler. Returns a positive integer alarm number
 * (can be used to cancel the pending timeout) or zero, if we are
 * returning from the timeout.
 **********************************************************************/
int
niceSetTimeout (unsigned int time)
{
  /* First, set the trap. This is the point the
   * niceNonLocalAlarmHandler () will unwind to (note:
   * niceNonLocalAlarmHandler () will always return FALSE) when the
   * alarm expires. If sigsetjmp instead returns 0, then it simply
   * means the trap is set. */
  if (sigsetjmp (timeoutenv, TRUE) == 0)
    return (niceSetAlarm (time, niceNonLocalAlarmHandler, FALSE, FALSE));

  /* The trap has sprung. Return FALSE = 0 to indicate that this is
   * the flow of execution after the timeout has expired. Note:
   * reversal of sigsetjmp return semantics to allow for duplexing the
   * alarm index with an indication that we are returning from the
   * longjmp. */
  return (FALSE);
}

/**********************************************************************
 * Reset a timeout. Assumes you had previously set and then cleared
 * the timeout and have therefore already set the appropriate non
 * local exit point.
 **********************************************************************/
int
niceResetTimeout (unsigned int time)
{
  return (niceSetAlarm (time, niceNonLocalAlarmHandler, FALSE, FALSE));
}

/**********************************************************************
 * Remove an alarm. Requires scanning down the freelist and snipping
 * out the appropriate cell. Be sure to increment the time field of
 * the next cell by the time field of the snipped cell (unless its the
 * first cell you are removing, in which case you add in the time
 * remaining to the alarm.
 **********************************************************************/
int
niceClearAlarm (int id)
{
#ifdef TMRDBG
  unsigned int tmp;
  AlarmCell *cell;
#endif
  /* How much time is left until the first pending alarm? */
  unsigned int timeLeft=alarm(0);

  /* Sanity check: no pending alarms and timeLeft is not zero. Someone
   * else is fooling around with your timer! */
  if (timeLeft == 0 && !alarms)
    return (ERROR);

  /* Update your rate estimate, if necessary. */
  if (++alarmOpCount > alarmOpInt)
    niceAlarmEstimate ();

  /* Remove the cell. */
  timeLeft = niceRemoveCell (id, timeLeft);

  /* Make estimated alarm rate correction. */
  if (timeLeft > 0 && !RANDINT(alarmOpRate))
    {
#if FALSE
      fprintf (stderr, "Alarm: adjusting.\n");
#endif      
      timeLeft--;
    }

#ifdef TMRDBG
  tmp = timeLeft;
  fprintf (stderr, "Alarm [-%d]", id);
  cell = alarms;
  while (cell)
    {
      fprintf (stderr, " [%d@%u:%u]%s%s", cell->id, tmp, cell->delay,
	       (cell->props & RECURRING)?"*":"", (cell->props & NONLOCAL)?"!":"");
      tmp = tmp + cell->delay;
      cell = cell->next;
    }
  fprintf (stderr, "\n");
#endif

  /* If there are pending alarms ready to fire, do so (TODO: what
   * happens if the alarm that fires is nonreturning? Is that a
   * problem? Or not?). Otherwise, reset the alarm (if necessary)
   * before you exit. If there are no remaining alarm signals, reset
   * alarmIndex so it doesn't overflow over time. */
  if (alarms && timeLeft == 0)
    niceFireAlarm (0);
  else if (alarms)
    alarm (timeLeft);
  else
    alarmIndex = 0;

  /* Return the id of the newly removed alarm. */
  return (id);
}

/**********************************************************************
 * niceAddCell places the given cell in the appropriate place on the
 * global alarms list, where the first cell on the alarms list (if
 * any) is scheduled to fire timeLeft seconds from now. Returns number
 * of seconds until first alarm on global alarm list.
 **********************************************************************/
unsigned int
niceAddCell (AlarmCell *new, unsigned int timeLeft)
{
  unsigned int now;
  AlarmCell *cell;

  if (!alarms || new->time < timeLeft)
    {
      /* No scheduled alarms or this alarm fires sooner than the
       * extant first alarm. Make the new cell the first cell. */
      new->next = alarms;
      alarms = new;
      new->delay = (new->next?(timeLeft - new->time):0);
      return (new->time);
    }

  /* OK, Alarms is non-NULL.  Find out where the new cell fits
   * in. When the loop finished, cell will point to the cell
   * immediately ahead of the first legal location for the new cell
   * (there may be multiple legal locations if there are multiple
   * cells firing at the same time), and now will tell you what time
   * cell is to fire. */
  cell = alarms;
  now = timeLeft;
  while (cell && cell->next && (now + cell->delay < new->time))
    {
      now = now + cell->delay;
      cell = cell->next;
    }
  /* At this point, cell points to the cell which should preceed the
   * new cell. Insert the new cell and make appropriate changes to the
   * time delays. */
  new->next = cell->next;
  cell->next = new;
  new->delay = (new->next?(now + cell->delay - new->time):0);
  cell->delay = new->time - now;

  /* One last adjustment: if the newly added cell has a nonreturning
   * handler, and the following cell is supposed to fire immediately,
   * we need to introduce an artificial delay so that the nonlocal
   * exit can be allowed to happen before firing the next alarm. */
  if (cell->props & NONLOCAL && cell->delay == 0)
    cell->delay = 1;

  /* Return amount of time to first cell firing, which in this case is
   * unchanged, since the first cell is unchanged. */
  return (timeLeft);
}

/**********************************************************************
 * niceRemoveCell removes the specified cell from the global alarms
 * list, where the first cell on the alarms list (if any) is scheduled
 * to fire timeLeft seconds from now. Returns number of seconds until
 * first alarm on alarms list, or ERROR if no alarms are left.
**********************************************************************/
int
niceRemoveCell (int id, int timeLeft)
{
  AlarmCell *cell = alarms, *last = NULL;

  /* Advance until cell points to the cell you want to remove. */
  while (cell && cell->id != id)
    {
      last = cell;
      cell = cell->next;
    }

  if (!cell)
    /* Cell not found. */
    return (timeLeft);
  else if (!last)
    {
      /* Deleting the first cell: increase timeLeft to reflect the
       * added wait time. */
      alarms = cell->next;
      timeLeft = timeLeft + cell->delay;
    }
  else if (cell && last)
    {
      /* Deleting either an internal cell or the last cell. If you've
       * just nuked the last cell, reset the new last cell's delay to
       * 0, just to be neat. */
      last->next = cell->next;
      last->delay = (cell->next)?(last->delay + cell->delay):0;
    }

  /* Recycle freed cell. */
  cell->next = freeAlarms;
  freeAlarms = cell;

  /* Return time left to first alarm. */
  return (timeLeft);
}

/**********************************************************************
 * This is the generic alarm handler used to handle every SIGALRM.
 * Dispatch on the alarm id to the appropriate alarm-specific handler,
 * and reset the alarm if necessary. Careful handling required if the
 * alarm-specific handler is nonreturning; also, an alarm-specific
 * handler may itself set an alarm, which requires the code be fully
 * reentrant.
 **********************************************************************/
static void
niceFireAlarm (int signum)
{
  unsigned int lead = 0;
  AlarmCell *cell;

#if FALSE
  fprintf (stderr, "Alarm [fire]:");
  cell = alarms;
  while (cell)
    {
      fprintf (stderr, " %s[%d@%u:%u]%s%s", 
	       (lead==0)?"+":"",
	       cell->id, lead, cell->delay,
	       (cell->props & RECURRING)?"*":"", (cell->props & NONLOCAL)?"!":"");
      lead = lead + cell->delay;
      cell = cell->next;
    }
  fprintf (stderr, "\n");
  lead = 0;
#endif

  /* Run all alarms that have expired: since you're here, it means the
   * first one has expired by definition. The only question is how
   * many of the leading cells are ready to go. Careful to reset cell
   * since the alarm list may reset while exec'ing the alarm
   * itself. */
  while ((cell = alarms) && (lead == 0))
    {
      /* Pop the cell and update lead time. */
      alarms = cell->next;
      lead = lead + cell->delay;
#ifdef TMRDBG
      fprintf (stderr, " [%d]%s%s", cell->id,
	       (cell->props & RECURRING)?"*":"", (cell->props & NONLOCAL)?"!":"");
#endif
      /* Reschedule or recycle the expired cell. */
      if (cell->props & RECURRING)
	{
	  /* Take this opportunity to reset the alarmIndex parameter
	   * if appropriate. Having a recurring alarm essentially
	   * means the alarm list will never be empty, which is your
	   * usual chance to reset the alarmIndex. So we need to be a
	   * bit more aggressive: if the alarm you are resetting is
	   * the only alarm, reset the alarmIndex to the recurring
	   * cell's ID, and start counting new alarms from there. */
	  if (!alarms)
	    alarmIndex = cell->id;	/* Keep current id for debugging ease. */
	  lead = niceAddCell (cell, lead);
	}
      else
	{
	  cell->next = freeAlarms;
	  freeAlarms = cell;
	}

      /* Restart the alarm before you call the handler. You need to do
       * this because the alarm handler may itself need to set an
       * alarm; furthermore, the handler may not return, in which case
       * you won't get a chance to set the alarm later. By resetting
       * the alarm, you ensure that any alarm manipulation that does
       * occur within the handler don't screw up. Of course, if there
       * are no other cells to fire it isn't an issue. 
       *
       * On tricky issue: if the current lead time is zero (that is,
       * the next cell is to fire right away) setting the alarm to 0
       * will not work correctly, because that is how we turn alarms
       * off. So we'll artificially delay the next cell by 1; since we
       * do adjust the timing of the alarms, the effect should even
       * out over time.
       *
       * If there are no alarm cells pending, reset alarmIndex (to
       * avoid an eventual identifier overflow). */
      if (alarms)
	alarm (MAX(lead, 1));
      else
	alarmIndex = 0;

      /* Call the handler. */
      (*(cell->handler)) ();

      /* OK, handler returned. Disable the alarm before you check for
       * additional cells to fire. */
      lead = alarm (0);
    }

  /* Fired all due alarms. Reset the alarm, or if there are no alarm
   * cells pending, reset alarmIndex (to avoid an eventual identifier
   * overflow). */
  if (alarms)
    alarm (MAX(lead,1));
  else
    alarmIndex = 0;
  return;
}

/**********************************************************************
 * niceAlarmEstimate is used to estimate the rate at which we perform
 * alarm operations, so that we can make corrections for too-frequent
 * ops.
 **********************************************************************/
void
niceAlarmEstimate ()
{
  time_t now = time(NULL);
  unsigned int elapsed;
  
  /* If this is the first time we've been called, there won't be any
   * cached time value yet. Cache the value, then exit: you'll get a
   * rate estimate next time. */
  if (alarmOpLastTime != 0)
    {
      /* Check the difference and change the alarm adjustment
       * parameters. alarmOpRate should always be greater than 1,
       * where 1 corresponds to "no time adjustment". */
      elapsed = MAX(1, (unsigned int) difftime (now, alarmOpLastTime));
      alarmOpRate = MAX(1, (alarmOpInt / elapsed));
  
      /* Adjust the counting interval so that it corresponds to about
       * 30 seconds at the current rate. This regulates how often the
       * alarm rate is calculated and adjusted. */
      alarmOpInt = 30*alarmOpRate;
#ifdef TMRDBG
      fprintf (stderr, "Alarm: rate=%u interval=%u\n", 
	       alarmOpRate, alarmOpInt);
#endif
    }
  /* Reset the counter and the time value. */
  alarmOpLastTime = now;
  alarmOpCount = 0;
  return;
}

/**********************************************************************
 * Default nonlocal timeout handler.  Uses nonlocal return; invoked
 * when blocked protocol is detected in order to achieve a clean
 * exit. Note that the jump environment is a single static global
 * variable (and thus private to this library), so I don't believe we
 * can handle nested timeouts with this default timeout alarm handler.
 * In other words, don't try setting a timeout inside
 * niceNonLocalAlarmHandler ()!
 **********************************************************************/
void
niceNonLocalAlarmHandler (int signum)
{
  /* All you need to do is invoke siglongjmp() with a non-zero (i.e.,
   * not FALSE) value. */
  siglongjmp (timeoutenv, TRUE);
}
#endif

/**********************************************************************
 * Message logging system. Based in part on a preliminary
 * implementation by James Hunsaker.
 **********************************************************************/
/* Tells you what level of message you are reporting. 0 level are
 * always reported (and cause termination). Otherwise, only errors of
 * the specified level or below are reported. */
int logLevel = NLOGFAT;

/* Tells you where each level message should be directed. */
int logDestination = NLOGSTDERR;

/**********************************************************************
 * Maps a message level, expressed as a collection of bits, to a
 * specified destination. Called at startup to specify where each type
 * of message should be directed. If this is not invoked, the default
 * is to issue only fatal messages on stderr.
 **********************************************************************/
void
niceLogSet (int level, int destination)
{
  logLevel = level;
  logDestination = destination;
  return;
}

/**********************************************************************
 * Issues given message to given destinations based on specified
 * message levels. The input levels is a bit pattern characterizing
 * the (possibly multiplex) severity of the message.
 **********************************************************************/
#ifndef WIN32
#define TMPBUFFERLEN 2048
#endif
void
niceLog (int severity, const char *format, ...)
{
  va_list argp;
#ifndef WIN32
  char temp[TMPBUFFERLEN+1];
#endif

  /* Output the message on each of the matching output streams, but
   * only if the indicated severity meets or exceeds (i.e., is less
   * than or equal to) the system setting. */
  if (logLevel >= severity)
    {
      if (logDestination & NLOGSTDOUT)
	{
	  va_start (argp, format);
	  vfprintf (stdout, format, argp);
	  va_end (argp);
	}
      if (logDestination & NLOGSTDERR)
	{
	  va_start (argp, format);
	  vfprintf (stderr, format, argp);
	  va_end (argp);
	}
#ifndef WIN32
      if (logDestination & NLOGSYSLOG)
	{
	  va_start (argp, format);
	  vsnprintf (temp, TMPBUFFERLEN, format, argp);
	  va_end (argp);
	  syslog (LOG_ERR, temp);
	}
#endif
      /* If fatal, exit; else simply return. */
      if (severity == NLOGFAT)
	exit (ERROR);
    }
  return;
}
