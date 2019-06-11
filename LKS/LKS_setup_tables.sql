DROP TABLE IF EXISTS Markers;
DROP TABLE IF EXISTS Maps;
DROP TABLE IF EXISTS MapMarkers;
DROP TABLE IF EXISTS Studies;
DROP TABLE IF EXISTS Analyses;
DROP TABLE IF EXISTS Pedigrees;
DROP TABLE IF EXISTS Positions;
DROP TABLE IF EXISTS PedigreePositions;
DROP TABLE IF EXISTS Regions;
DROP TABLE IF EXISTS Servers;
DROP VIEW IF EXISTS servers;
DROP TABLE IF EXISTS DModelParts;
DROP TABLE IF EXISTS QModelParts;
DROP TABLE IF EXISTS TP2MP;
DROP TABLE IF EXISTS Models;
DROP TABLE IF EXISTS LGModels;
DROP TABLE IF EXISTS MarkerSetLikelihood;
DROP TABLE IF EXISTS MarkerSetLikelihood_MCMC;
DROP TABLE IF EXISTS ServerPedigrees;


CREATE TABLE Markers (
StudyId int NOT NULL,
Name varchar(16) NOT NULL COMMENT 'Frequencies are not standard, but rather specific to pedigree map',
ChromosomeNo int NULL,
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
PRIMARY KEY (StudyId, Name)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Maps (
MapId int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
MapScale char(1) NOT NULL COMMENT 'Will be (K)osambi or (H)aldane',
ReferenceFlag tinyint NOT NULL COMMENT '1-ReferenceMap 0-StudyMap',
Description varchar(128) COMMENT 'Use and preface with a MEANINGFUL file name',
-- Really can't use this since StudyId won't exist initially... CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
CONSTRAINT UNIQUE KEY MinimalMaps (StudyId, ReferenceFlag, MapScale, Description),
PRIMARY KEY (MapId)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE MapMarkers (
StudyId int NOT NULL,
MapId int NOT NULL,
MarkerName varchar(16) NOT NULL,
RefPosCM real NULL COMMENT 'Start-scaling and trait-applicable range reference position',
AvePosCM real NULL COMMENT 'Start-scaling marker genotype position',
NextRefPosCM real NULL COMMENT 'End of trait-applicable range for this marker',
ScaleRefPosCM real NULL COMMENT 'End-scaling marker reference position',
ScaleAvePosCM real NULL COMMENT 'End-scaling marker genotype position',
Scale real NULL,
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
-- CONSTRAINT FOREIGN KEY (MapId) references Maps (MapId),
-- CONSTRAINT FOREIGN KEY (MarkerName) references Markers (Name),
PRIMARY KEY (StudyId, MapId, MarkerName)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Studies (
StudyId int NOT NULL AUTO_INCREMENT,
StudyLabel varchar(64) NOT NULL COMMENT 'Should be concisely descriptive and uniq',
ReferenceMapId int NOT NULL DEFAULT 0 COMMENT 'Merge results onto this map',
LiabilityClassCnt int NOT NULL COMMENT 'Reasonably 1, 2 or 3',
ImprintingFlag char(1) NOT NULL COMMENT 'Y/N',
Description varchar(128),
PendingWorkFlag char(1) DEFAULT NULL COMMENT 'Null is indeterminate status, Y = work available, N = no work available, D = study done',
InsertTime timestamp DEFAULT CURRENT_TIMESTAMP,
-- CONSTRAINT FOREIGN KEY (ReferenceMapId) references Maps (MapId),
CONSTRAINT UNIQUE KEY (StudyLabel), 
PRIMARY KEY (StudyId)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Analyses (
StudyId int(11) NOT NULL,
Uniquey varbinary(767),
PedigreeRegEx varchar(1024) NOT NULL,
PedigreeNotRegEx varchar(1024) NOT NULL DEFAULT 'XYZZY',
AnalysisId int(11) NOT NULL AUTO_INCREMENT,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
PRIMARY KEY (AnalysisId),
UNIQUE KEY (StudyId, Uniquey)
) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Pedigrees (
StudyId int NOT NULL,
PedigreeSId varchar(16) NOT NULL COMMENT 'String ID from file',
GenotypeMapId int NULL COMMENT 'Convert from this map to ReferenceMapId',
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
-- CONSTRAINT FOREIGN KEY (GenotypeMapId) references Maps (MapId),
PRIMARY KEY (StudyId, PedigreeSId)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Positions (
StudyId int NOT NULL,
ChromosomeNo int NOT NULL,
RefTraitPosCM real NOT NULL,
PPL double precision NULL,
PositionStatus char(1) NULL,
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
PRIMARY KEY (StudyId, ChromosomeNo, RefTraitPosCM)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE PedigreePositions (
PedPosId int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
PedigreeSId varchar(16) NOT NULL,
ChromosomeNo int NOT NULL,
RefTraitPosCM real NOT NULL,
PedTraitPosCM real NULL,
MarkerCount int DEFAULT 0 COMMENT 'Minimum desired MP marker count',
FreeModels int DEFAULT 0 COMMENT 'Amount of available work for this ped/pos',
PendingLikelihoods int DEFAULT 0 COMMENT 'Amount of incomplete work for this ped/pos',
SingleModelEstimate int NULL COMMENT 'Estimated non-polynomial runtime of models at this ped/pos (from SMRT run)',
SingleModelRuntime int NULL COMMENT 'Actual peak non-polynomial runtime of models at this ped/pos (initially an estimate)',
PRIMARY KEY (PedPosId),
-- CONSTRAINT FOREIGN KEY (StudyId, PedigreeSId) references Pedigrees (StudyId, PedigreeSId),
INDEX (StudyId, PedigreeSId),
-- CONSTRAINT FOREIGN KEY (StudyId, ChromosomeNo, RefTraitPosCM) references Positions (StudyId, ChromosomeNo, RefTraitPosCM),
INDEX (StudyId, ChromosomeNo, PedTraitPosCM),
INDEX (StudyId, ChromosomeNo, PedTraitPosCM, MarkerCount, FreeModels),
INDEX (StudyId, ChromosomeNo, RefTraitPosCM)
) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Regions (
RegionId int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
AnalysisId int NOT NULL,
ChromosomeNo int NOT NULL,
RefTraitPosCM real NOT NULL,
RegionNo int NOT NULL,
ParentRegionNo int NULL COMMENT 'Split-from RegionNo (or 0)',
ParentRegionError double NULL COMMENT 'Error that triggered the split',
ParentRegionSplitDir int NULL COMMENT 'For Sang to say how parent split',
InsertTime timestamp DEFAULT CURRENT_TIMESTAMP,
CONSTRAINT UNIQUE KEY MinimalRegions (StudyId, AnalysisId, ChromosomeNo, RefTraitPosCM, RegionNo),
PRIMARY KEY (RegionId),
-- CONSTRAINT FOREIGN KEY (StudyId, ChromosomeNo, RefTraitPosCM) references Positions (StudyId, ChromosomeNo, RefTraitPosCM),
INDEX (StudyId, AnalysisId, ChromosomeNo, RefTraitPosCM, RegionNo)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Servers (
ServerId int NOT NULL AUTO_INCREMENT,
ConnectionId bigint(21) unsigned NOT NULL COMMENT 'Generated by MySQL',
HostName varchar(32) NOT NULL,
ProcessId int NOT NULL,
ListenerSocketId int NULL COMMENT 'For server groups with a centralized work disperser',
KeepAliveFlag int NOT NULL COMMENT 'Server will start with 1, set to 0 to shutdown',
StudyId int NOT NULL,
PedigreeRegEx varchar(1024) NOT NULL,
PedigreeNotRegEx varchar(1024) NOT NULL DEFAULT 'XYZZY',
ChromosomeNo int NOT NULL,
Algorithm varchar(2) NOT NULL COMMENT 'LG or ES',
MarkerCount int NOT NULL,
ProgramVersion varchar(32) NULL,
StartTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
LastHeartbeat timestamp NULL COMMENT 'Currently updated with each cache refill',
CurrentPedPosId int NULL COMMENT 'PedPosId we are working on (for combined likelihood only)',
CurrentLimit int NULL COMMENT 'Max number of models per pass for this PedPosId',
StopTime timestamp NULL,
ExitStatus int NULL,

-- MCMC
  SampleIdStart int NULL COMMENT 'Applicable to MCMC likelihood servers only', 
  SampleIdEnd int NULL COMMENT 'Applicable to MCMC likelihood servers only',
 
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
INDEX (StudyId),
PRIMARY KEY (ServerId)) ENGINE=InnoDB comment='$Id$';

CREATE VIEW servers as select
ServerId, ConnectionId, HostName, ProcessId, StudyId, StartTime, LastHeartbeat, CurrentPedPosId, CurrentLimit, StopTime, ExitStatus from Servers;

CREATE TABLE DModelParts (
MPId int NOT NULL AUTO_INCREMENT,
DGF decimal(32,30) NOT NULL,
BigPen decimal(32,30) NOT NULL,
BigLittlePen decimal(32,30) NOT NULL,
LittleBigPen decimal(32,30) NOT NULL,
LittlePen decimal(32,30) NOT NULL,
PRIMARY KEY (MPId),
UNIQUE KEY ModelByValues (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE QModelParts (
MPId int NOT NULL AUTO_INCREMENT,
DGF decimal(32,30) NOT NULL,
BigMean decimal(32,30) NOT NULL,
BigLittleMean decimal(32,30) NOT NULL,
LittleBigMean decimal(32,30) NOT NULL,
LittleMean decimal(32,30) NOT NULL,
BigSD decimal(32,30) NOT NULL,
BigLittleSD decimal(32,30) NOT NULL,
LittleBigSD decimal(32,30) NOT NULL,
LittleSD decimal(32,30) NOT NULL,
Threshold decimal(32,30) NOT NULL,
PRIMARY KEY (MPId),
UNIQUE KEY ModelByValues (DGF, BigMean, BigLittleMean, LittleBigMean, LittleMean,
	BigSD, BigLittleSD, LittleBigSD, LittleSD, Threshold)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE TP2MP (
ModelId int NOT NULL,
WeightedLRComponent double DEFAULT NULL,
RuntimeCostSec int NULL,
PRIMARY KEY (ModelId)) Engine=InnoDB comment='$Id$';

CREATE TABLE Models (
ModelId int NOT NULL AUTO_INCREMENT,
PedPosId int NOT NULL,
LC1MPId int NOT NULL,
LC2MPId int NOT NULL,
LC3MPId int NOT NULL,
MarkerCount int DEFAULT 0,
ServerId int NULL,
StartTime timestamp NULL,
Likelihood double NULL,
RuntimeCostSec int NULL,
EndTime timestamp NULL,
-- CONSTRAINT FOREIGN KEY (PedPosId) references PedigreePositions (PedPosId),
-- INDEX (PedPosId),
-- CONSTRAINT FOREIGN KEY (LC1MPId) references DModelParts (MPId),
INDEX (LC1MPId),
-- CONSTRAINT FOREIGN KEY (LC2MPId) references DModelParts (MPId),
-- CONSTRAINT FOREIGN KEY (LC3MPId) references DModelParts (MPId),
-- CONSTRAINT FOREIGN KEY (ServerId) references Servers (ServerId),
INDEX (ServerId),
UNIQUE KEY (ModelId),
-- INDEX (PedPosId, LC1MPId, LC2MPId, LC3MPId, ServerId),
PRIMARY KEY (PedPosID, LC1MPID, LC2MPId, LC3MPId, MarkerCount)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE LGModels (
LGModelID int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
PedigreeRegEx varchar(1024) NOT NULL,
PedigreeNotRegEx varchar(1024) NOT NULL DEFAULT 'XYZZY',
ServerId int,
LC1MPId int,
LC2MPId int,
LC3MPId int,
PRIMARY KEY (LGModelId)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE MarkerSetLikelihood (
  MarkerSetId int NOT NULL AUTO_INCREMENT,
  PedPosId int NOT NULL, 
  MarkerCount int NULL,
  ServerId int NULL, 
  StartTime timestamp NULL, 
  Likelihood double NULL,
  RuntimeCostSec int NULL, 
  EndTime timestamp NULL,
-- CONSTRAINT FOREIGN KEY (PedPosId) references PedigreePositions (PedPosId),
INDEX (PedPosId),
-- CONSTRAINT FOREIGN KEY (ServerId) references Servers (ServerId),
INDEX (ServerId),
INDEX(PedPosID, MarkerCount),
PRIMARY KEY (MarkerSetId)
) ENGINE=InnoDB comment='$Id$';

CREATE TABLE MarkerSetLikelihood_MCMC (
  MarkerSetId int NOT NULL, 
  SampleId int NOT NULL,
  Likelihood double NULL,
  INDEX(MarkerSetId), 
  PRIMARY KEY (MarkerSetId, SampleId)
) ENGINE=InnoDB comment='$Id$';

CREATE TABLE ServerPedigrees (
  ServerId int NOT NULL, 
  PedigreeSId varchar(16) NOT NULL,

  PRIMARY KEY (ServerId, PedigreeSId)
) ENGINE=InnoDB comment='$Id$';


-- These are a couple of views used for testing purposes; test cases can pull
-- these into files and then do a comparison to make sure we haven't broken
-- anything. Only two so far, but hoping to add more.
DROP VIEW IF EXISTS models_likelihood_list;
DROP VIEW IF EXISTS mcmc_likelihood_list;

-- Models likelihood list - a simple join on the Models table that shows pedpos
-- info for each model.
-- FIXME: this doesn't really come up with unique rows with non-Likelihood
-- columns, which makes it potentially problematic for comparisons. The Right
-- Way to do that would be to also do joins with DModelParts and possibly also
-- QModelParts, but I'm not sure how to usefully represent that. -jo
CREATE VIEW models_likelihood_list AS
    SELECT PP.PedigreeSId,
            PP.ChromosomeNo,
            PP.RefTraitPosCM,
            M.Likelihood
    FROM Models M
    LEFT JOIN PedigreePositions PP ON M.PedPosId = PP.PedPosId
    ORDER BY PedigreeSId, ChromosomeNo, RefTraitPosCM, Likelihood;

-- MCMC likelihood list - shows per-sample per-pedigree-position likelihoods.
CREATE VIEW mcmc_likelihood_list AS
    SELECT PP.PedigreeSId,
            PP.ChromosomeNo,
            PP.RefTraitPosCM,
            MSLMCMC.SampleId,
            MSLMCMC.Likelihood
    FROM MarkerSetLikelihood_MCMC MSLMCMC
    LEFT JOIN MarkerSetLikelihood MSL
            ON MSLMCMC.MarkerSetId = MSL.MarkerSetId
    LEFT JOIN PedigreePositions PP ON MSL.PedPosId = PP.PedPosId
    ORDER BY PedigreeSId, ChromosomeNo, RefTraitPosCM, SampleId;


-- "Unused" tables follow
-- 
-- Technically, these tables aren't actually regularly used at the moment. They
-- are retained nonetheless because while some team members here work best with
-- comments in code, other team members here work best with the presence of
-- tables in the database. So we leave them in the schema, and change whether
-- or not they're getting used in LKS_setup_trigger_proc.sql.
DROP TABLE IF EXISTS RegionModels;
DROP TABLE IF EXISTS HundredBlockSizingRuns;
DROP TABLE IF EXISTS SingleSizingRuns;
DROP TABLE IF EXISTS SingleModelRuntimes;
DROP TABLE IF EXISTS Diag;

-- Indicates which Regions are affected by a particular ModelId
CREATE TABLE RegionModels (
RegionId int NOT NULL,
ModelId int NOT NULL,
PRIMARY KEY (RegionId, ModelId),
INDEX (ModelId)
) ENGINE=InnoDB comment='$Id$' AUTO_INCREMENT=1110;

CREATE TABLE HundredBlockSizingRuns (
Directory varchar(255) NOT NULL,
HundredBlock int(11) NOT NULL AUTO_INCREMENT,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
PRIMARY KEY (Directory),
KEY `HundredBlock` (HundredBlock)
) ENGINE=InnoDB AUTO_INCREMENT=1110 comment='$Id$';

CREATE TABLE SingleSizingRuns (
Directory varchar(255) NOT NULL,
StudyId int(11) NOT NULL AUTO_INCREMENT,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
PRIMARY KEY (Directory),
KEY `StudyId` (StudyId)
) ENGINE=InnoDB AUTO_INCREMENT=100000 comment='$Id$';

-- Must be independent of PedigreePositions both for multiple MarkerCount issues and map-merge issues
CREATE TABLE SingleModelRuntimes (
StudyId int NOT NULL,
PedigreeSId varchar(16) NOT NULL,
ChromosomeNo int NOT NULL,
PedTraitPosCM real NOT NULL,
MarkerCount int DEFAULT 0,
SingleModelRuntime int DEFAULT NULL,
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
PRIMARY KEY (StudyId, PedigreeSId, ChromosomeNo, PedTraitPosCM, MarkerCount)) ENGINE=InnoDB comment='$Id$';

CREATE TABLE Diag (
DiagId int NOT NULL AUTO_INCREMENT,
ConnectionId bigint(21) unsigned,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
Message text,
PRIMARY KEY (DiagId)) ENGINE=InnoDB comment='$Id$';
