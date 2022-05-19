-- Copyright (C) 2009, 2022 Mathematical Medicine LLC
-- 
-- This program is free software: you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the Free
-- Software Foundation, either version 3 of the License, or (at your option)
-- any later version.
-- 
-- You should have received a copy of the GNU General Public License along
-- with this program. If not, see <https://www.gnu.org/licenses/>.
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



