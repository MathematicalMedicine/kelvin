Create table RunLog (
entryNo integer auto_increment primary key,
fromNode varchar(32),
logEntry varchar(255),
updateDate timestamp default current_timestamp on update current_timestamp);

Create table RunLogProcessing (lastEntryNo integer);

Create table RunStatus (
nodeName varchar(32),
PID integer,
startTime datetime,
kelvinVersion varchar(64),
configFile varchar(128),
runType varchar(128),
finishStats varchar(128),
endTime datetime,
primary key (nodeName, PID, startTime));



