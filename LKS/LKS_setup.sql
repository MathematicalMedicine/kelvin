-- Do things like: Insert into Diag (Message) values (Concat('Leaving GetDLikelihood with Likelihood of ', @outLikelihood)); for out-of-band diagnostics..

DROP TRIGGER IF EXISTS IncModels;
DROP TRIGGER IF EXISTS DecModels;
DROP TABLE IF EXISTS TP2MP;
DROP TABLE IF EXISTS Models;
DROP TABLE IF EXISTS LGModels;
DROP TABLE IF EXISTS QModelParts;
DROP TABLE IF EXISTS DModelParts;
DROP TABLE IF EXISTS Servers;
DROP VIEW IF EXISTS servers;
DROP TABLE IF EXISTS Regions;
DROP TABLE IF EXISTS SingleModelRuntimes;
DROP TABLE IF EXISTS PedigreePositions;
DROP TABLE IF EXISTS Positions;
DROP TABLE IF EXISTS Pedigrees;
DROP TABLE IF EXISTS Studies;
DROP TABLE IF EXISTS S2R;
DROP TABLE IF EXISTS MapMarkers;
DROP TABLE IF EXISTS Maps;
DROP TABLE IF EXISTS Markers;
DROP TABLE IF EXISTS Diag;
DROP TABLE IF EXISTS SingleSizingRuns;
DROP TABLE IF EXISTS HundredBlockSizingRuns;

CREATE TABLE HundredBlockSizingRuns (
Directory varchar(255) NOT NULL,
HundredBlock int(11) NOT NULL AUTO_INCREMENT,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
PRIMARY KEY (Directory),
KEY `HundredBlock` (HundredBlock)
) ENGINE=InnoDB AUTO_INCREMENT=1110;

CREATE TABLE SingleSizingRuns (
Directory varchar(255) NOT NULL,
StudyId int(11) NOT NULL AUTO_INCREMENT,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
PRIMARY KEY (Directory),
KEY `StudyId` (StudyId)
) ENGINE=InnoDB AUTO_INCREMENT=100000;

CREATE TABLE Diag (
DiagId int NOT NULL AUTO_INCREMENT,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
Message varchar(120) NULL,
PRIMARY KEY (DiagId)) ENGINE=InnoDB;

CREATE TABLE Markers (
StudyId int NOT NULL,
Name varchar(16) NOT NULL COMMENT 'Frequencies are not standard, but rather specific to pedigree map',
ChromosomeNo int NULL,
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
PRIMARY KEY (StudyId, Name)) ENGINE=InnoDB;

CREATE TABLE Maps (
MapId int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
MapScale char(1) NOT NULL COMMENT 'Will be (K)osambi or (H)aldane',
Description varchar(128) COMMENT 'Use and preface with a MEANINGFUL file name',
-- Really can't use this since StudyId won't exist initially... CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
CONSTRAINT UNIQUE KEY MinimalMaps (StudyId, MapScale, Description),
PRIMARY KEY (MapId)) ENGINE=InnoDB;

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
PRIMARY KEY (StudyId, MapId, MarkerName)) ENGINE=InnoDB;

CREATE TABLE Studies (
StudyId int NOT NULL AUTO_INCREMENT,
ReferenceMapId int NOT NULL COMMENT 'Merge results onto this map',
LiabilityClassCnt int NOT NULL COMMENT 'Reasonably 1, 2 or 3',
ImprintingFlag char(1) NOT NULL COMMENT 'Y/N',
Description varchar(128),
PendingWorkFlag char(1) DEFAULT NULL COMMENT 'Null is indeterminate status, Y = work available, N = no work available, D = study done',
-- CONSTRAINT FOREIGN KEY (ReferenceMapId) references Maps (MapId),
PRIMARY KEY (StudyId)) ENGINE=InnoDB;

CREATE TABLE Pedigrees (
StudyId int NOT NULL,
PedigreeSId varchar(16) NOT NULL COMMENT 'String ID from file',
GenotypeMapId int NULL COMMENT 'Convert from this map to ReferenceMapId',
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
-- CONSTRAINT FOREIGN KEY (GenotypeMapId) references Maps (MapId),
PRIMARY KEY (StudyId, PedigreeSId)) ENGINE=InnoDB;

CREATE TABLE Positions (
StudyId int NOT NULL,
ChromosomeNo int NOT NULL,
RefTraitPosCM real NOT NULL,
PPL double precision NULL,
PositionStatus char(1) NULL,
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
PRIMARY KEY (StudyId, ChromosomeNo, RefTraitPosCM)) ENGINE=InnoDB;

CREATE TABLE PedigreePositions (
PedPosId int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
PedigreeSId varchar(16) NOT NULL,
ChromosomeNo int NOT NULL,
RefTraitPosCM real NOT NULL,
PedTraitPosCM real NULL,
MarkerCount int NULL COMMENT 'Minimum desired MP marker count',
FreeModels int DEFAULT 0 COMMENT 'Amount of available work for this ped/pos',
SingleModelEstimate int NULL COMMENT 'Estimated non-polynomial runtime of models at this ped/pos (from SMRT run)',
SingleModelRuntime int NULL COMMENT 'Actual peak non-polynomial runtime of models at this ped/pos (initially an estimate)',
PRIMARY KEY (PedPosId),
-- CONSTRAINT FOREIGN KEY (StudyId, PedigreeSId) references Pedigrees (StudyId, PedigreeSId),
INDEX (StudyId, PedigreeSId),
-- CONSTRAINT FOREIGN KEY (StudyId, ChromosomeNo, RefTraitPosCM) references Positions (StudyId, ChromosomeNo, RefTraitPosCM),
INDEX (StudyId, ChromosomeNo, RefTraitPosCM)) ENGINE=InnoDB;

-- Must be independent of PedigreePositions both for multiple MarkerCount issues and map-merge issues
CREATE TABLE SingleModelRuntimes (
StudyId int DEFAULT NULL,
PedigreeSId varchar(16) NOT NULL,
ChromosomeNo int NOT NULL,
PedTraitPosCM real NOT NULL,
MarkerCount int DEFAULT NULL,
SingleModelRuntime int DEFAULT NULL,
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
PRIMARY KEY (StudyId, PedigreeSId, ChromosomeNo, PedTraitPosCM, MarkerCount)) ENGINE=InnoDB;

CREATE TABLE Regions (
RegionId int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
ChromosomeNo int NOT NULL,
RefTraitPosCM real NOT NULL,
RegionNo int NOT NULL,
ParentRegionNo int NULL COMMENT 'Split-from RegionNo (or 0)',
ParentRegionError double NULL COMMENT 'Error that triggered the split',
ParentRegionSplitDir int NULL COMMENT 'For Sang to say how parent split',
PendingLikelihoods int DEFAULT 0 COMMENT 'How many remain to be computed',
CONSTRAINT UNIQUE KEY MinimalRegions (StudyId, ChromosomeNo, RefTraitPosCM, RegionNo),
PRIMARY KEY (RegionId),
-- CONSTRAINT FOREIGN KEY (StudyId, ChromosomeNo, RefTraitPosCM) references Positions (StudyId, ChromosomeNo, RefTraitPosCM),
INDEX (StudyId, ChromosomeNo, RefTraitPosCM, RegionNo)) ENGINE=InnoDB;

CREATE TABLE Servers (
ServerId int NOT NULL AUTO_INCREMENT,
ConnectionId int NOT NULL COMMENT 'Generated by MySQL',
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
-- CONSTRAINT FOREIGN KEY (StudyId) references Studies (StudyId),
INDEX (StudyId),
PRIMARY KEY (ServerId)) ENGINE=InnoDB;

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
UNIQUE KEY ModelByValues (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)) ENGINE=InnoDB;

CREATE TABLE QModelParts (
MPId int NOT NULL AUTO_INCREMENT,
DGF decimal(32,30) NOT NULL,
BigPenMean decimal(32,30) NOT NULL,
BigLittlePenMean decimal(32,30) NOT NULL,
LittleBigPenMean decimal(32,30) NOT NULL,
LittlePenMean decimal(32,30) NOT NULL,
BigPenSD decimal(32,30) NOT NULL,
BigLittlePenSD decimal(32,30) NOT NULL,
LittleBigPenSD decimal(32,30) NOT NULL,
LittlePenSD decimal(32,30) NOT NULL,
Threshold decimal(32,30) NOT NULL,
PRIMARY KEY (MPId),
UNIQUE KEY ModelByValues (DGF, BigPenMean, BigLittlePenMean, LittleBigPenMean, LittlePenMean,
	BigPenSD, BigLittlePenSD, LittleBigPenSD, LittlePenSD, Threshold)) ENGINE=InnoDB;

CREATE TABLE TP2MP (
ModelId int NOT NULL,
WeightedLRComponent double DEFAULT NULL,
RuntimeCostSec int NULL,
PRIMARY KEY (ModelId)) Engine=InnoDB;

CREATE TABLE Models (
ModelId int NOT NULL AUTO_INCREMENT,
PedPosId int NOT NULL,
LC1MPId int NOT NULL,
LC2MPId int NOT NULL,
LC3MPId int NOT NULL,
MarkerCount int NULL,
RegionId int NOT NULL,
ServerId int NULL,
StartTime timestamp NULL,
Likelihood double NULL,
RuntimeCostSec int NULL,
EndTime timestamp NULL,
-- CONSTRAINT FOREIGN KEY (PedPosId) references PedigreePositions (PedPosId),
INDEX (PedPosId),
-- CONSTRAINT FOREIGN KEY (LC1MPId) references DModelParts (MPId),
INDEX (LC1MPId),
-- CONSTRAINT FOREIGN KEY (LC2MPId) references DModelParts (MPId),
-- CONSTRAINT FOREIGN KEY (LC3MPId) references DModelParts (MPId),
-- CONSTRAINT FOREIGN KEY (RegionId) references Regions (RegionId),
INDEX (RegionId),
-- CONSTRAINT FOREIGN KEY (ServerId) references Servers (ServerId),
INDEX (ServerId),
UNIQUE KEY (ModelId),
PRIMARY KEY (PedPosID, LC1MPID, LC2MPId, LC3MPId, MarkerCount)) ENGINE=InnoDB;

CREATE TABLE LGModels (
LGModelID int NOT NULL AUTO_INCREMENT,
StudyId int NOT NULL,
ServerId int,
LC1MPId int,
LC2MPId int,
LC3MPId int,
PRIMARY KEY (LGModelId)) ENGINE=InnoDB;

Delimiter //
CREATE TRIGGER DecModels AFTER UPDATE ON Models FOR EACH ROW
BEGIN
  Update PedigreePositions set FreeModels = FreeModels - 1 where
    PedPosId = NEW.PedPosId AND
    OLD.ServerId IS NULL AND
    NEW.ServerId IS NOT NULL;
  Update Regions set PendingLikelihoods = PendingLikelihoods - 1 where
    RegionId = NEW.RegionId AND
    OLD.Likelihood IS NULL AND
    NEW.Likelihood IS NOT NULL;
  Update PedigreePositions set FreeModels = FreeModels + 1 where
    PedPosId = NEW.PedPosId AND
    OLD.ServerId IS NOT NULL AND
    NEW.ServerId IS NULL;
  Update Regions set PendingLikelihoods = PendingLikelihoods + 1 where
    RegionId = NEW.RegionId AND
    OLD.Likelihood IS NOT NULL AND
    NEW.Likelihood IS NULL;
END;
//
CREATE TRIGGER IncModels AFTER INSERT ON Models FOR EACH ROW
BEGIN
  Update PedigreePositions set FreeModels = FreeModels + 1 where
    PedPosId = NEW.PedPosId AND
    NEW.ServerId IS NULL AND
    NEW.Likelihood IS NULL;
  Update Regions set PendingLikelihoods = PendingLikelihoods + 1 where
    RegionId = NEW.RegionId AND
    NEW.Likelihood IS NULL;
END;
//
Delimiter ;

DROP PROCEDURE IF EXISTS BadScaling;
DELIMITER //

-- BadScaling is called after the addition of new trait positions or new pedigree genotype
-- maps. It (re)populates the MapMarkers scaling information and PedigreePositions PedTraitPosCM
-- column. It uses a bad linear interpolation approach for marker-bounded positions, and a
-- smarter extrapolation approach for the rest.

CREATE PROCEDURE BadScaling (IN inStudyId int)
BEGIN
  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';
  DECLARE no_rows_indicator INT DEFAULT 0;
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

WholeThing: LOOP

  -- Add a new temporary reference map for every pedigree map with a MapId of pedigree MapId+100000
  -- This allows us to ensure that only markers common to both reference and pedigree maps are
  -- present.
  Insert ignore into MapMarkers (StudyId, MapId, MarkerName, RefPosCM, AvePosCM, NextRefPosCM)
    Select a.StudyId, b.MapId+100000, a.MarkerName, a.RefPosCM, a.AvePosCM, a.NextRefPosCM
    from MapMarkers a, MapMarkers b, Studies c where
    c.StudyId = inStudyId AND
    a.StudyId = c.StudyId AND b.StudyId = c.StudyId AND
    a.MapId = c.ReferenceMapId AND
    b.MapId <> c.ReferenceMapId AND b.MapId < 100000 AND
    a.MarkerName = b.MarkerName;

  -- Make sure the reference map has a RefPosCM value
  Update MapMarkers a, Studies b set a.RefPosCM = a.AvePosCM where a.MapId = b.ReferenceMapId;

  -- Get (personal) reference map positions on every pedigree map.
  Update MapMarkers a, Studies b, MapMarkers c set a.RefPosCM = c.AvePosCM where
	a.StudyId = inStudyId AND a.StudyId = b.StudyId AND b.StudyId = c.StudyId AND
	c.MapId = a.MapId+100000 AND a.MapId < 100000 AND
	a.MarkerName = c.MarkerName;

  -- Insert a new temporary leftmost marker to permit scaling of positions to the left of the first marker
  Insert ignore into MapMarkers (StudyId, MapId, MarkerName, RefPosCM, AvePosCM, NextRefPosCM)
    Select StudyId, MapId, 'Dummy', min(RefPosCM)-320, min(RefPosCM)-320, min(RefPosCM)
      from MapMarkers where StudyId = inStudyId group by MapId;

  -- Scaling table has next reference position...
  Create temporary table Scaling
	Select a.StudyId, a.MapId, a.MarkerName, min(b.RefPosCM) 'NextPosCM' from 
		MapMarkers a, MapMarkers b where 
		b.MapId = a.MapId+100000 AND
		a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
		a.RefPosCM < b.RefPosCM AND a.MarkerName <> b.MarkerName 
		group by a.MapId, a.MarkerName;

  Update MapMarkers a, Scaling b
	set a.NextRefPosCM = b.NextPosCM, a.ScaleRefPosCM = b.NextPosCM where
	a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
	a.MapId = b.MapId AND a.MarkerName = b.MarkerName;

  Drop temporary table Scaling;

  -- Scaling table has next pedigree position...
  Create temporary table Scaling
	Select a.StudyId, a.MapId, a.MarkerName, min(b.AvePosCM) 'NextPosCM' from 
		MapMarkers a, MapMarkers b where
		a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
		a.AvePosCM < b.AvePosCM AND
		a.MarkerName <> b.MarkerName AND a.MapId = b.MapId
		group by a.MapId, a.MarkerName;

  Update MapMarkers a, Scaling b
	set a.ScaleAvePosCM = b.NextPosCM where
		a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
		a.MapId = b. MapId AND a.MarkerName = b.MarkerName;

  Update MapMarkers set Scale = (ScaleAvePosCM - AvePosCM) / (ScaleRefPosCM - RefPosCM) where
	StudyId = inStudyId;

  Drop temporary table Scaling;

  Downscale: LOOP

    SET no_rows_indicator = 0;
    Select MapId, MarkerName into @MapId, @MarkerName from MapMarkers where
	StudyId = inStudyId AND
	(Scale > 10.0 OR Scale < 0.1) AND Scale IS NOT NULL limit 1;
    IF no_rows_indicator THEN
      LEAVE Downscale;
    END IF;

    Create temporary table Scaling
	Select a.StudyId, a.MapId, a.MarkerName, min(b.ScaleRefPosCM) 'NextPosCM' from
	        MapMarkers a, MapMarkers b where
		a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
		(a.Scale > 10.0 OR a.Scale < 0.1) AND a.Scale IS NOT NULL AND
		a.ScaleRefPosCM < b.ScaleRefPosCM AND a.MarkerName <> b.MarkerName AND a.MapId = b.MapId
	        group by a.MapId, a.MarkerName;
    Update MapMarkers a, Scaling b set a.ScaleRefPosCM = b.NextPosCM where
	a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
	a.MapId = b. MapId AND a.MarkerName = b.MarkerName;
    Drop temporary table Scaling;

    Create temporary table Scaling
	Select a.StudyId, a.MapId, a.MarkerName, min(b.ScaleAvePosCM) 'NextPosCM' from
	        MapMarkers a, MapMarkers b where
		a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
		(a.Scale > 10.0 OR a.Scale < 0.1) AND a.Scale IS NOT NULL AND
		a.ScaleAvePosCM < b.ScaleAvePosCM AND a.MarkerName <> b.MarkerName AND a.MapId = b.MapId
	        group by a.MapId, a.MarkerName;
    Update MapMarkers a, Scaling b set a.ScaleAvePosCM = b.NextPosCM where
	a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
	a.MapId = b. MapId AND a.MarkerName = b.MarkerName;
    Update MapMarkers set Scale = (ScaleAvePosCM - AvePosCM) / (ScaleRefPosCM - RefPosCM) where
	StudyId = inStudyId;
    Drop temporary table Scaling;

  END LOOP Downscale;

  -- Deal with the rightmost marker by setting scale to the next rightmost.
  Update MapMarkers a, MapMarkers b set a.Scale = b.Scale where
	a.StudyId = inStudyId AND a.StudyId = b.StudyId AND
	a.Scale IS NULL AND a.MapId = b.MapId AND a.RefPosCM = b.NextRefPosCM;

  -- Indicate that the scaling associated with the remaining (rightmost) markers is good clear
  -- out to beyond the end of any chromosome.
  Update MapMarkers set NextRefPosCM = 999 where
	StudyId = inStudyId AND
	NextRefPosCM is NULL;

  -- Scale all of the positions using the factors associated with the left marker.x
  Update PedigreePositions a, Pedigrees b, MapMarkers c 
	set a.PedTraitPosCM = c.AvePosCM + ((a.RefTraitPosCM - c.RefPosCM) * c.Scale) where 
	a.StudyId = inStudyId AND a.StudyId = b.StudyId AND b.StudyId = c.StudyId AND
	a.PedigreeSId = b.PedigreeSId AND 
	b.GenotypeMapId = c.MapId AND 
	a.RefTraitPosCM >= c.RefPosCM AND a.RefTraitPosCM < c.NextRefPosCM;

  Delete from MapMarkers where MarkerName = 'Dummy';

  -- Yeah, this is just a diagnostic, a way of bailing early to test incrementally.
  Leave WholeThing;

END LOOP WholeThing;

END;
//
Delimiter ;

DROP PROCEDURE IF EXISTS GetDLikelihood;
DELIMITER //

-- GetLikelihood is called by the kelvin client to request and retrieve Likelihoods. All possible fixed-model
-- parameters are passed in the call. They should be -1 if they're N/A. If a corresponding Likelihood
-- is found, it is returned, otherwise the call is treated as an asynchronous request and a NULL
-- is returned for the Likelihood. Generally speaking, actual results are only guarenteed once the
-- pertinent region is marked as complete in the Regions table.

CREATE PROCEDURE GetDLikelihood (
 IN inPedPosId int, IN inDGF real, 
 IN inLC1BigPen real, IN inLC1BigLittlePen real, IN inLC1LittleBigPen real, IN inLC1LittlePen real,
 IN inLC2BigPen real, IN inLC2BigLittlePen real, IN inLC2LittleBigPen real, IN inLC2LittlePen real,
 IN inLC3BigPen real, IN inLC3BigLittlePen real, IN inLC3LittleBigPen real, IN inLC3LittlePen real,
 IN inRegionNo int, IN inParentRegionNo int, IN inParentRegionError real, IN inParentRegionSplitDir int,
 OUT outRegionId int, OUT outMarkerCount int, OUT outLikelihood real
)
BEGIN
    DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';
    DECLARE WNE_indicator INT DEFAULT 0; -- For when we know the Likelihood _will_ not exist

    DECLARE no_rows_indicator INT DEFAULT 0;
    DECLARE no_rows CONDITION FOR 1329;
    DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

    Start transaction;

    -- There's always at least one liability class...
    SET no_rows_indicator = 0;
    Select MPId into @LC1MPId from DModelParts where
    DGF = convert(inDGF, DECIMAL(32,30)) AND BigPen = convert(inLC1BigPen, DECIMAL(32,30)) AND BigLittlePen = convert(inLC1BigLittlePen, DECIMAL(32,30)) AND 
    LittleBigPen = convert(inLC1LittleBigPen, DECIMAL(32,30)) AND LittlePen = convert(inLC1LittlePen, DECIMAL(32,30));

    IF no_rows_indicator THEN
      Insert into DModelParts (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)
      values (inDGF, inLC1BigPen, inLC1BigLittlePen, inLC1LittleBigPen, inLC1LittlePen);
      Select LAST_INSERT_ID() INTO @LC1MPId;
      SET WNE_indicator = 1;
    END IF;

    SET no_rows_indicator = 0;
    Select MPId into @LC2MPId from DModelParts where
    DGF = convert(inDGF, DECIMAL(32,30)) AND BigPen = convert(inLC2BigPen, DECIMAL(32,30)) AND BigLittlePen = convert(inLC2BigLittlePen, DECIMAL(32,30)) AND 
    LittleBigPen = convert(inLC2LittleBigPen, DECIMAL(32,30)) AND LittlePen = convert(inLC2LittlePen, DECIMAL(32,30));
    
    IF no_rows_indicator THEN
      Insert into DModelParts (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)
      values (inDGF, inLC2BigPen, inLC2BigLittlePen, inLC2LittleBigPen, inLC2LittlePen);
      Select LAST_INSERT_ID() INTO @LC2MPId;
      SET WNE_indicator = 1;
    END IF;

    SET no_rows_indicator = 0;
    Select MPId into @LC3MPId from DModelParts where
    DGF = convert(inDGF, DECIMAL(32,30)) AND BigPen = convert(inLC3BigPen, DECIMAL(32,30)) AND BigLittlePen = convert(inLC3BigLittlePen, DECIMAL(32,30)) AND 
    LittleBigPen = convert(inLC3LittleBigPen, DECIMAL(32,30)) AND LittlePen = convert(inLC3LittlePen, DECIMAL(32,30));
    
    IF no_rows_indicator THEN
      Insert into DModelParts (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)
      values (inDGF, inLC3BigPen, inLC3BigLittlePen, inLC3LittleBigPen, inLC3LittlePen);
      Select LAST_INSERT_ID() INTO @LC3MPId;
      SET WNE_indicator = 1;
    END IF;

    -- Now we've got all 3 MPIds, get the RegionId and find or add the row in Models

    Select StudyId, ChromosomeNo, RefTraitPosCM INTO @StudyId, @ChromosomeNo, @RefTraitPosCM
    from PedigreePositions where PedPosId = inPedPosId;

    Set no_rows_indicator = 0;
    Select RegionId into outRegionId from Regions where
    StudyId = @StudyId AND ChromosomeNo = @ChromosomeNo AND 
    RefTraitPosCM = @RefTraitPosCM AND RegionNo = inRegionNo;

    If no_rows_indicator THEN
      Insert into Regions (StudyId, ChromosomeNo, RefTraitPosCM, RegionNo, ParentRegionNo, ParentRegionError, ParentRegionSplitDir) values
      (@StudyId, @ChromosomeNo, @RefTraitPosCM, inRegionNo, inParentRegionNo, inParentRegionError, inParentRegionSplitDir);
      Select LAST_INSERT_ID() INTO outRegionId;
    END IF;

    IF WNE_indicator THEN
      Insert into Models (PedPosId, LC1MPId, LC2MPId, LC3MPId, RegionId) values
	(inPedPosId, @LC1MPId, @LC2MPId, @LC3MPId, outRegionId);
    ELSE
      -- Go for the Likelihood! If it doesn't exist, request it and return the NULL.
      SET no_rows_indicator = 0;
      Select Likelihood, RegionId, MarkerCount into outLikelihood, outRegionId, outMarkerCount from
	Models where PedPosId = inPedPosId AND LC1MPId = @LC1MPId AND
	  LC2MPId = @LC2MPId AND LC3MPId = @LC3MPId order by MarkerCount desc limit 1;
	
      IF no_rows_indicator THEN
        Insert into Models (PedPosId, LC1MPId, LC2MPId, LC3MPId, RegionId) values
	  (inPedPosId, @LC1MPId, @LC2MPId, @LC3MPId, outRegionId);
      END IF;
    END IF;

    Commit;
END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetDParts;
DELIMITER //

CREATE PROCEDURE GetDParts (
  IN inLC1MPId int, inLC2MPId int, inLC3MPId int,
  OUT outDGF real,
  OUT outLC1BP real, OUT outLC1BLP real, OUT outLC1LBP real, OUT outLC1LP real,
  OUT outLC2BP real, OUT outLC2BLP real, OUT outLC2LBP real, OUT outLC2LP real,
  OUT outLC3BP real, OUT outLC3BLP real, OUT outLC3LBP real, OUT outLC3LP real
)
BEGIN
  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

  -- Better to return the results as a result set than out parameters (which suck!)
  Select LC1.DGF,
	LC1.BigPen, LC1.BigLittlePen, LC1.LittleBigPen, LC1.LittlePen,
	LC2.BigPen, LC2.BigLittlePen, LC2.LittleBigPen, LC2.LittlePen,
	LC3.BigPen, LC3.BigLittlePen, LC3.LittleBigPen, LC3.LittlePen
  into
	outDGF,
	outLC1BP, outLC1BLP, outLC1LBP, outLC1LP, 
	outLC2BP, outLC2BLP, outLC2LBP, outLC2LP, 
	outLC3BP, outLC3BLP, outLC3LBP, outLC3LP
  from DModelParts LC1, DModelParts LC2, DModelParts LC3
  where LC1.MPId = inLC1MPId AND LC2.MPId = inLC2MPId AND LC3.MPId = inLC3MPId;
    
--  Insert into Diag (Message) values (Concat('GetDParts: returning ', convert(outDGF,char),', ',
--	convert(outLC1BP,char),', ', convert(outLC1BLP,char),', ', convert(outLC1LBP,char),', ', convert(outLC1LP,char),', ', 
--	convert(outLC2BP,char),', ', convert(outLC2BLP,char),', ', convert(outLC2LBP,char),', ', convert(outLC2LP,char),', ', 
--	convert(outLC3BP,char),', ', convert(outLC3BLP,char),', ', convert(outLC3LBP,char),', ', convert(outLC3LP,char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetWork;
DELIMITER //

-- GetWork is called by the AltL servers to acquire an available fixed model that matches their
-- pedigree/position pattern. Only QT/DT-independent information is provided. OUT parameters
-- need not be selected as they will be returned in the only result set.

CREATE PROCEDURE GetWork (
 IN inServerId int, IN inLowPosition real, IN inHighPosition real, IN inLocusListType int,
 OUT outPedPosId int, OUT outPedigreeSId varchar(32), OUT outPedTraitPosCM real,
 OUT outLC1MPId int, OUT outLC2MPId int, OUT outLC3MPId int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';
  DECLARE no_rows_indicator INT DEFAULT 0;
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

--  Insert into Diag (Message) values (Concat('GetWork: called w/ ', convert(inServerId,char),', ',convert(inLowPosition,char),
--	', ',convert(inHighPosition,char)));

  Create temporary table if not exists CachedWork (
	WorkId int auto_increment,
	PedPosId int, PedigreeSId varchar(32), PedTraitPosCM real,
	LC1MPId int, LC2MPId int, LC3MPId int,
	ServerId int, MarkerCount int, RegionId int,
	PRIMARY KEY (WorkId));

-- See if our cached work table has a row for us to return...

  SET no_rows_indicator = 0;
  Select WorkId, PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId
  into @WorkId, outPedPosId, outPedigreeSId, outPedTraitPosCM, outLC1MPId, outLC2MPId, outLC3MPId
  from CachedWork order by PedigreeSId limit 1;

  IF no_rows_indicator THEN

    -- Someday we may want to verify that we didn't accidently change positions here and orphan a lot of work
  
    IF inLocusListType = 1 THEN
      call CacheMarkerWork(inServerId, inLowPosition, inHighPosition, @ResultRows);
    ELSE
      IF inLocusListType = 2 THEN
        call CacheTraitWork(inServerId, @ResultRows);
      ELSE
        IF inLocusListType = 3 THEN
          call CacheCombinedWork(inServerId, inLowPosition, inHighPosition, @ResultRows);
        ELSE
          IF inLocusListType = 100 THEN
            call Cache2ptInitialWork(inServerId, inLowPosition, inHighPosition, @ResultRows);
          ELSE
            IF inLocusListType = 101 THEN
              call Cache2ptFinalWork(inServerId, inLowPosition, inHighPosition, @ResultRows);
            ELSE
	      Insert into Diag (Message) values (Concat('GetWork: unexpected inLocusListType of ',convert(inLocusListType,char)));
            END IF;
          END IF;
        END IF;
      END IF;
    END IF;
 
    -- Now try again...

    SET no_rows_indicator = 0;
    Select WorkId, PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId
    into @WorkId, outPedPosId, outPedigreeSId, outPedTraitPosCM, outLC1MPId, outLC2MPId, outLC3MPId
    from CachedWork order by PedigreeSId limit 1;

  END IF;

  IF no_rows_indicator THEN
--    Insert into Diag (Message) values ('GetWork: no work found at all');
    Truncate table CachedWork;
  ELSE
--    Insert into Diag (Message) values (Concat('GetWork: returning ', convert(outPedPosId,char),', ', outPedigreeSId, ', ',
--    convert(outPedTraitPosCM,char),', ', convert(outLC1MPId,char),', ', convert(outLC2MPId,char),', ', convert(outLC3MPId,char)));
    Delete from CachedWork where WorkId = @WorkId;
  END IF;

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS Cache2ptInitialWork;
DELIMITER //

-- Cache2ptInitialWork is like the old GetWork but returns multiple rows in a temporary table (specific to connection) named CachedWork.

CREATE PROCEDURE Cache2ptInitialWork (
 IN inServerId int, IN inLowPosition real, IN inHighPosition real, OUT outResultRows int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

  Insert into Diag (Message) values (Concat('Cache2ptInitialWork: called w/ ', convert(inServerId,char),', ',convert(inLowPosition,char),
	', ',convert(inHighPosition,char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  Cache2ptInitialWork: LOOP

    Select KeepAliveFlag into @KeepAliveFlag from Servers where ServerId = inServerId;
    IF @KeepAliveFlag = 0 THEN
      LEAVE Cache2ptInitialWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select fresh (no ServerId or StartTime), work which is generally to the right of the current marker
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount, RegionId)
      select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount, 0
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId IS NULL
      limit 50 for update; -- No ordering for best performance
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN
	Insert into Diag (Message) values ('Cache2ptInitialWork: no work found at all');
      LEAVE Cache2ptInitialWork;
    ELSE
      -- Work we found has not been done at all, update the rows
      Insert into Diag (Message) values (Concat('Cache2ptInitialWork: found ', convert(outResultRows,char), ' rows of undone work, marking'));
      Select count(*) from CachedWork into outResultRows;
      Update Models a, CachedWork b  set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId IS NULL;
      LEAVE Cache2ptInitialWork;
    END IF;

  END LOOP Cache2ptInitialWork;

  Commit;

  Insert into Diag (Message) values (Concat('Cache2ptInitialWork: returning ', convert(outResultRows,char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS Cache2ptFinalWork;
DELIMITER //

-- Cache2ptFinalWork is like the old GetWork but returns multiple rows in a temporary table (specific to connection) named CachedWork.

CREATE PROCEDURE Cache2ptFinalWork (
 IN inServerId int, IN inLowPosition real, IN inHighPosition real, OUT outResultRows int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

  Insert into Diag (Message) values (Concat('Cache2ptFinalWork: called w/ ', convert(inServerId,char),', ',convert(inLowPosition,char),
	', ',convert(inHighPosition,char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  Cache2ptFinalWork: LOOP

    Select KeepAliveFlag into @KeepAliveFlag from Servers where ServerId = inServerId;
    IF @KeepAliveFlag = 0 THEN
      LEAVE Cache2ptFinalWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select started but incomplete (has our ServerId, probably an entry in TP2MP) generally to the left of the current marker
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount, RegionId)
      select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount, 0
      from
	PedigreePositions PP, Models M, Servers S, TP2MP T
      where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId = S.ServerId AND
	M.ModelId = T.ModelId AND
	M.EndTime IS NULL
      limit 50 for update; -- No ordering for best performance
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    LEAVE Cache2ptFinalWork;

  END LOOP Cache2ptFinalWork;

  Commit;

  Insert into Diag (Message) values (Concat('Cache2ptFinalWork: returning ', convert(outResultRows,char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS CacheCombinedWork;
DELIMITER //

-- CacheCombinedWork is like the old GetWork but returns multiple rows in a temporary table (specific to connection) named CachedWork.

CREATE PROCEDURE CacheCombinedWork (
 IN inServerId int, IN inLowPosition real, IN inHighPosition real, OUT outResultRows int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

  -- Insert into Diag (Message) values (Concat('CacheCombinedWork: called w/ ', convert(inServerId,char),', ',convert(inLowPosition,char),
  --	', ',convert(inHighPosition,char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;

  Truncate table CachedWork;

  Create temporary table if not exists ExcludedPedPosIds (PedPosId int);

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite, so be careful!
  CacheCombinedWork: LOOP

    Select KeepAliveFlag into @KeepAliveFlag from Servers where ServerId = inServerId;
    IF @KeepAliveFlag = 0 THEN
      Insert into Diag (Message) values ('CacheCombinedWork: no longer kept alive');
      LEAVE CacheCombinedWork;
    END IF;

    SET outResultRows = 0;
    SET @candidatePedPosId = NULL;
    SET @candidateLimit = 51;

    Start transaction;

    -- Pick a random candidate PedigreePosition and limit Models to that PedPosId so we get a tolerable total runtime.
    -- Note that a server will have done all Trait and Marker Set models FIRST for any position, so we don't have to
    -- exclude them from our work selection.

    Select PP.PedPosId, PP.SingleModelRuntime, PP.FreeModels
      into @candidatePedPosId, @candidateSMRT, @candidateFreeModels
      from PedigreePositions PP, Servers S where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.PedPosId NOT IN (Select PedPosId from ExcludedPedPosIds) AND
	PP.FreeModels > 0 AND
	PP.MarkerCount = S.MarkerCount AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition
        order by RAND() limit 1 for update;
    	
    IF @candidatePedPosId IS NULL THEN
      -- Insert into Diag (Message) values (Concat('CacheCombinedWork: no work found for PedPosId ', convert(@candidatePedPosId,char)));
      SET outResultRows = 0;
      LEAVE CacheCombinedWork;
    --    ELSE
    --  Insert into Diag (Message) values (Concat('CacheCombinedWork: work with SMRT of ', convert(@candidateSMRT,char), ' found for PedPosId ', convert(@candidatePedPosId,char)));
    END IF;

    IF @candidateSMRT IS NULL THEN
      SET @candidateLimit = 5;
    ELSE
      IF @candidateSMRT = 0 THEN
        SET @candidateLimit = 51;
      ELSE
        SET @candidateLimit = convert((60*60) / @candidateSMRT, decimal(5,0));
        IF @candidateLimit < 1 THEN
          SET @candidateSMRT = 0;
          SET @candidateLimit = 1;
        END IF;
        IF @candidateLimit > 51 THEN
          SET @candidateLimit = 51;
        END IF;
      END IF;
    END IF;
    Insert into Diag (Message) values (Concat('CacheCombinedWork: work for PedPosId ', convert(@candidatePedPosId,char), ' with SMRT ', convert(@candidateSMRT,char), ' limited to ', convert(@candidateLimit,char)));
    Update Servers set CurrentPedPosId = @candidatePedPosId, CurrentLimit = @candidateLimit  where ServerId = inServerId;

    SET @dSString = Concat(
    	'Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount, RegionId) ',
	'select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount, 0 ',
	'from PedigreePositions PP, Models M, Servers S ',
	'where S.ServerId = ', convert(inServerId,char),
	' AND PP.PedPosId = ', convert(@candidatePedPosId,char),
	' AND M.PedPosId = ', convert(@candidatePedPosId,char),
	' AND M.ServerId IS NULL limit ',
	convert(@candidateLimit,char), ' for update;');
    PREPARE dSHandle from @dSString;
    EXECUTE dSHandle;
    DEALLOCATE PREPARE dSHandle;

    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN
      Insert into ExcludedPedPosIds (PedPosId) value (@candidatePedPosId);
      Select count(*) from ExcludedPedPosIds into @foo;
      Insert into Diag (Message) values (Concat('CacheCombinedWork: no actual work (not ', convert(@candidateFreeModels,char), ') found for ServerId ', convert(inServerId,char), ', PedPosId ', convert(@candidatePedPosId,char), ', excluded ', convert(@foo,char), ', looping'));
      -- We don't leave, we try again with a different PedPosId (hopefully)
      Commit;
    ELSE
      -- Insert into Diag (Message) values (Concat('CacheCombinedWork: found ', convert(outResultRows,char), ' rows of undone work, marking'));
      Update Models a, CachedWork b  set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId IS NULL;
      LEAVE CacheCombinedWork;
    END IF;

  END LOOP CacheCombinedWork;

  Drop table ExcludedPedPosIds;
  Commit;

  -- Insert into Diag (Message) values (Concat('CacheCombinedWork: returning ', convert(outResultRows,char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS CacheMarkerWork;
DELIMITER //

-- CacheMarkerWork is like the old GetWork but returns multiple rows in a temporary table (specific to connection) named CachedWork.

CREATE PROCEDURE CacheMarkerWork (
 IN inServerId int, IN inLowPosition real, IN inHighPosition real, OUT outResultRows int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

--  Insert into Diag (Message) values (Concat('CacheMarkerWork: called w/ ', convert(inServerId,char),', ',convert(inLowPosition,char),
--	', ',convert(inHighPosition,char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  CacheMarkerWork: LOOP

    Select KeepAliveFlag into @KeepAliveFlag from Servers where ServerId = inServerId;
    IF @KeepAliveFlag = 0 THEN
      LEAVE CacheMarkerWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select first undone (no ServerId), and then if need be, less quality (lesser MarkerCount) work.
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount, RegionId)
      select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount, 0
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.LC1MPId in (Select distinct MPId from DModelParts where DGF = -1) AND
	M.ServerId IS NULL
      limit 50 for update; -- No ordering for best performance
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN

--      Insert into Diag (Message) values ('CacheMarkerWork: found no undone work, checking lesser quality work');

      SET outResultRows = 0;
      Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount, RegionId)
        select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M1.LC1MPId, M1.LC2MPId, M1.LC3MPId, M1.ServerId, S.MarkerCount, M1.RegionId
        from
	  PedigreePositions PP, Models M1, Servers S
        where
	  S.ServerId = inServerId AND
	  PP.StudyId = S.StudyId AND
	  PP.MarkerCount <= S.MarkerCount AND
	  PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	  PP.ChromosomeNo = S.ChromosomeNo AND
	  PP.PedTraitPosCM > inLowPosition AND
	  PP.PedTraitPosCM <= inHighPosition AND
	  PP.PedPosId = M1.PedPosId AND
	  M1.LC1MPId in (Select distinct MPId from DModelParts where DGF = -1) AND
	  M1.MarkerCount < S.MarkerCount AND
	  NOT EXISTS (
		Select * from Models M2 where 
		M1.PedPosId = M2.PedPosId AND
		M1.LC1MPId = M2.LC1MPId AND M1.LC2MPId = M2.LC2MPId AND M1.LC3MPId = M2.LC3MPId AND
		M2.MarkerCount = S.MarkerCount
		)
	  limit 50 for update;
      Select count(*) from CachedWork into outResultRows;

      IF outResultRows = 0 THEN
--	Insert into Diag (Message) values ('CacheMarkerWork: no work found at all');
        LEAVE CacheMarkerWork;
      ELSE
--	Insert into Diag (Message) values (Concat('CacheMarkerWork: found ', convert(outResultRows,char), ' rows of lesser-quality, prepping by inserting new Models row'));
        Insert into Models (PedPosId, LC1MPId, LC2MPId, LC3MPId, RegionId, ServerId, MarkerCount, StartTime)
	 select PedPosId, LC1MPId, LC2MPId, LC3MPId, RegionId, ServerId, MarkerCount, CURRENT_TIMESTAMP from CachedWork;
        LEAVE CacheMarkerWork;
      END IF;

    ELSE
      -- Work we found has not been done at all, update the rows
--      Insert into Diag (Message) values (Concat('CacheMarkerWork: found ', convert(outResultRows,char), ' rows of undone work, marking'));
      Select count(*) from CachedWork into outResultRows;
      Update Models a, CachedWork b  set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId IS NULL;
      LEAVE CacheMarkerWork;
    END IF;

  END LOOP CacheMarkerWork;

  Commit;

--  Insert into Diag (Message) values (Concat('CacheMarkerWork: returning ', convert(outResultRows,char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS CacheTraitWork;
DELIMITER //

-- CacheTraitWork is like the old GetWork but returns multiple rows in a temporary table (specific to connection) named CachedWork.

CREATE PROCEDURE CacheTraitWork (
 IN inServerId int, OUT outResultRows int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

--  Insert into Diag (Message) values (Concat('CacheTraitWork: called w/ ', convert(inServerId,char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  CacheTraitWork: LOOP

    Select KeepAliveFlag into @KeepAliveFlag from Servers where ServerId = inServerId;
    IF @KeepAliveFlag = 0 THEN
      LEAVE CacheTraitWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select first undone (no ServerId), and then if need be, less quality (lesser MarkerCount) work.
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount, RegionId)
      select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount, 0
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.RefTraitPosCM = -9999.99 AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId IS NULL
      limit 50 for update; -- No ordering for best performance
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN
--      Insert into Diag (Message) values ('CacheTraitWork: found no undone work (none at all for traits!)');
      LEAVE CacheTraitWork;
    ELSE
      -- Work we found has not been done at all, update the rows
--      Insert into Diag (Message) values (Concat('CacheTraitWork: found ', convert(outResultRows,char), ' rows of undone work, marking'));
      Select count(*) from CachedWork into outResultRows;
      Update Models a, CachedWork b  set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId IS NULL;
      LEAVE CacheTraitWork;
    END IF;

  END LOOP CacheTraitWork;

  Commit;

  --  Insert into Diag (Message) values (Concat('CacheTraitWork: returning ', convert(outResultRows,char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS CountWork;
DELIMITER //

CREATE PROCEDURE CountWork (
 IN inServerId int, IN inLowPosition real, IN inHighPosition real,
 OUT outWorkCount int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';
  DECLARE no_rows_indicator INT DEFAULT 0;
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

--  Insert into Diag (Message) values (Concat('CountWork: called w/ ', convert(inServerId,char),', ',convert(inLowPosition,char),
--	', ',convert(inHighPosition,char)));

  -- Select first undone (no ServerId), and then if need be, less quality (lesser MarkerCount) work.

  SET no_rows_indicator = 0;
  Select count(*) into outWorkCount
  from
      PedigreePositions PP, Models M, Servers S
  where
	S.ServerId = inServerId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId IS NULL;

--  Insert into Diag (Message) values (Concat('CountWork: found ', convert(outWorkCount,char),' undone work'));
  IF no_rows_indicator THEN
    Select count(*) into outWorkCount
    from
        PedigreePositions PP, Models M1, Servers S
    where
        S.ServerId = 1 AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND
	PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M1.PedPosId AND
	M1.MarkerCount < S.MarkerCount AND
	NOT EXISTS (
		Select * from Models M2 where 
		M1.PedPosId = M2.PedPosId AND
		M1.LC1MPId = M2.LC1MPId AND M1.LC2MPId = M2.LC2MPId AND M1.LC3MPId = M2.LC3MPId AND
		M2.MarkerCount = S.MarkerCount
		);
--    Insert into Diag (Message) values (Concat('CountWork: found ', convert(outWorkCount,char),' lower-quality  work'));
  END IF;

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS PutWork;
DELIMITER //

-- PutWork is called by Likelihood servers to put the calculated Likelihood score in the Models table.

CREATE PROCEDURE PutWork (
  IN inServerId int, IN inPedPosId int, IN inLC1MPId int, IN inLC2MPId int, IN inLC3MPId int, IN inMarkerCount int, IN inLikelihood double, IN inRuntimeCostSec int
)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: Likelihood_server.sql 118 2010-10-07 15:24:05Z whv001 $';

--  Insert into Diag (Message) values (Concat('PutWork: called w/ ', convert(inServerId,char), ', ', convert(inPedPosId,char),', ',
--	convert(inLC1MPId,char),', ', convert(inLC2MPId,char),', ', convert(inLC3MPId,char), ', ',
--	convert(inMarkerCount,char),', ', convert(inLikelihood,char),', ', convert(inRuntimeCostSec,char)));

  Set @newSMRT = 0;

  IF inMarkerCount = 100 THEN
    -- 2pt, only have the initial first half, hold onto it.
    Select ModelId INTO @outModelId from Models where
	PedPosId = inPedPosId AND
	LC1MPId = inLC1MPId AND
	LC2MPId = inLC2MPId AND
	LC3MPId = inLC3MPId AND
	ServerId = inServerId;
    Insert into TP2MP (ModelId, WeightedLRComponent, RuntimeCostSec) values (@outModelId, inLikelihood, inRuntimeCostSec);
  ELSE
    IF inMarkerCount = 101 THEN
      -- 2pt, got the final half, combine them and truely finish the job.
      Select ModelId INTO @outModelId from Models where
	PedPosId = inPedPosId AND
	LC1MPId = inLC1MPId AND
	LC2MPId = inLC2MPId AND
	LC3MPId = inLC3MPId AND
	ServerId = inServerId;
      Select WeightedLRComponent, RuntimeCostSec into @outWeightedLRComponent, @newSMRT from TP2MP where ModelId = @outModelId;
      Insert into Diag (Message) values (Concat('PutWork: Adding inLikelihood ', convert(inLikelihood,char), ' to existing ', convert(@outWeightedLRComponent,char)));
      Update Models set Likelihood = inLikelihood+@outWeightedLRComponent, RuntimeCostSec = inRuntimeCostSec+@newSMRT,
	MarkerCount = inMarkerCount, EndTime = CURRENT_TIMESTAMP where ModelId = @outModelId;
    ELSE
      -- Classic multipoint
      Update Models set Likelihood = inLikelihood, RuntimeCostSec = inRuntimeCostSec, MarkerCount = inMarkerCount, EndTime = CURRENT_TIMESTAMP where
	PedPosId = inPedPosId AND
	LC1MPId = inLC1MPId AND
	LC2MPId = inLC2MPId AND
	LC3MPId = inLC3MPId AND
	ServerId = inServerId;
    END IF;
  END IF;

  -- Now the hard part. If this model took an unexpectedly long time, we need to update the PedPosId with that fact, and
  -- we want to jettison the remaining models because  we don't want to be holding lots of work.
  -- If you wonder why we have to do this, look at the variation on SingleModelRuntimes using a query like this:
  -- Select PedPosId, count(*), max(RuntimeCostSec), min(RuntimeCostSec), avg(RuntimeCostSec), stddev_pop(RuntimeCostSec) from Models where PedPosId in (Select PedPosId from PedigreePositions where StudyId = 8) group by PedPosId order by max(RuntimeCostSec) desc limit 20;
  IF inRuntimeCostSec > 60 THEN
    Select SingleModelRuntime into @oldSMRT from PedigreePositions where PedPosId = inPedPosId;
    IF @oldSMRT IS NULL THEN
      SET @oldSMRT = 0;
    END IF;
    IF (inRuntimeCostSec * 0.70) > @oldSMRT THEN
      -- Set the SingleModelRuntime to the new value
      Update PedigreePositions set SingleModelRuntime = inRuntimeCostSec where PedPosId = inPedPosId;
      Insert into Diag (Message) values (Concat('PutWork: Reset SMRT due to ', convert(inRuntimeCostSec,char), 's model for PedPosId ', convert(inPedPosId,char), ', which had ', convert(@oldSMRT,char), 's SMRT'));
    END IF;
    -- We can't just look at the oldSMRT because we could be holding work for which some other server just updated the cost.
    Select inRuntimeCostSec * count(*) into @oldLimit from CachedWork;
    IF @oldLimit > (60*60*0.70) THEN
      -- Release the models we haven't processed
      Update Models a, CachedWork b  set a.ServerId = NULL, a.MarkerCount = b.MarkerCount, a.StartTime = NULL where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId = inServerId;
      -- Clear-out everything so we get a new, fresh bunch of PROPERLY SIZED work
      Delete from CachedWork;
      Insert into Diag (Message) values (Concat('PutWork: flushed remaining work due to ', convert(inRuntimeCostSec,char), 's model for PedPosId ', convert(inPedPosId,char), ', which had ', convert(@oldSMRT,char), 's SMRT'));
    END IF;
  END IF;

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS SetDummyNullLikelihood;
DELIMITER //

CREATE PROCEDURE SetDummyNullLikelihood (IN inServerId int)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

-- Set all trait and marker likelhood Models for pedigrees handled by this ServerId to be finished and have an Likelihood of 1.0

  DECLARE inModelId int;
  DECLARE EmptyCursor BOOLEAN DEFAULT 0;
  DECLARE DummyMarkerModels CURSOR FOR
      Select M.ModelId
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND PP.StudyId = S.StudyId AND PP.FreeModels > 0 AND PP.MarkerCount = 1 AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedPosId = M.PedPosId AND M.LC1MPId in (Select distinct MPId from DModelParts where DGF = -1) AND M.ServerId IS NULL;
  DECLARE DummyTraitModels CURSOR FOR
      Select M.ModelId
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND PP.StudyId = S.StudyId AND PP.FreeModels > 0 AND PP.MarkerCount = 1 AND
	PP.PedigreeSId RLIKE S.PedigreeRegEx AND PP.PedigreeSId NOT RLIKE S.PedigreeNotRegEx AND PP.ChromosomeNo = S.ChromosomeNo AND
	PP.RefTraitPosCM = -9999.99 AND PP.PedPosId = M.PedPosId AND M.ServerId IS NULL;
  DECLARE CONTINUE HANDLER FOR SQLSTATE '02000' SET EmptyCursor = 1;

  Insert into Diag (Message) values (Concat('SetDummyNullLikelihood: called w/ ', convert(inServerId,char)));

  OPEN DummyMarkerModels;
  REPEAT
    FETCH DummyMarkerModels into inModelId;
    Update Models set ServerId = inServerId, StartTime = CURRENT_TIMESTAMP, Likelihood = 1.0, EndTime = CURRENT_TIMESTAMP where ModelId = inModelId;
  UNTIL EmptyCursor END REPEAT;
  CLOSE DummyMarkerModels;

  -- Then the trait Models

  SET EmptyCursor = 0;
  OPEN DummyTraitModels;
  REPEAT
    FETCH DummyTraitModels into inModelId;
    Update Models set ServerId = inServerId, StartTime = CURRENT_TIMESTAMP, Likelihood = 1.0, EndTime = CURRENT_TIMESTAMP where ModelId = inModelId;
  UNTIL EmptyCursor END REPEAT;
  CLOSE DummyTraitModels;

END;
//
Delimiter ;

DROP PROCEDURE IF EXISTS CleanOrphans;
DELIMITER //

CREATE PROCEDURE CleanOrphans (IN inStudyId int)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

-- Put orphaned Models back up for adoption.

  DECLARE inModelId int;
  DECLARE inPedPosId int;
  DECLARE orphanCount int DEFAULT 0;

  DECLARE EmptyCursor BOOLEAN DEFAULT 0;
  DECLARE OrphanModels CURSOR FOR
    Select ModelId, PedPosId from Models where
	StartTime is NOT NULL AND Likelihood is NULL AND
	PedPosId in (Select PedPosId from PedigreePositions where StudyId = inStudyId) AND
	ServerId NOT IN (Select ServerId from Servers where ExitStatus IS NULL);
  DECLARE CONTINUE HANDLER FOR SQLSTATE '02000' SET EmptyCursor = 1;

  OPEN OrphanModels;

  FETCH OrphanModels into inModelId, inPedPosId;
  WHILE EmptyCursor <> 1 DO

    UPDATE Models set ServerId = NULL, StartTime = NULL where ModelId = inModelId;
    UPDATE PedigreePositions set FreeModels = FreeModels + 1 where PedPosId = inPedPosId;
    DELETE from TP2MP where ModelId = inModelId;
    SET orphanCount = orphanCount+1;
    FETCH OrphanModels into inModelId, inPedPosId;

  END WHILE;

  CLOSE OrphanModels;

  Select orphanCount 'Orphans Freed';

END;
//
Delimiter ;

DROP PROCEDURE IF EXISTS CleanStudy;
DELIMITER //

CREATE PROCEDURE CleanStudy (IN inStudyId int)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

-- Remove all rows associated with a study

  Delete from TP2MP where ModelId in (Select ModelId from Models where PedPosId in (Select PedPosId from PedigreePositions where StudyId = inStudyId));
  Delete from Models where PedPosId in (Select PedPosId from PedigreePositions where StudyId = inStudyId);
  Delete from Servers where StudyId = inStudyId;
  Delete from Regions where StudyId = inStudyId;
  Delete from PedigreePositions where StudyId = inStudyId;
  Delete from Positions where StudyId = inStudyId;
  Delete from Pedigrees where StudyId = inStudyId;
  Delete from MapMarkers where StudyId = inStudyId;
  Delete from Markers where StudyId = inStudyId;
  Delete from Maps where StudyId = inStudyId;
  Delete from Studies where StudyId = inStudyId;

END;
//
Delimiter ;

DROP PROCEDURE IF EXISTS MoveStudy;
DELIMITER //

CREATE PROCEDURE MoveStudy (IN inFromStudyId int, IN inToStudyId int)
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

-- Move all rows associated with a StudyId to a new StudyId without violating integrity constraints

  -- Studies row could already exist
  Update ignore Studies set StudyId = inToStudyId where StudyId = inFromStudyId;
  Delete ignore from Studies where StudyId = inFromStudyId;
  Update ignore Maps set StudyId = inToStudyId where StudyId = inFromStudyId;
  Delete ignore from Maps where StudyId = inFromStudyId;
  Update ignore Markers set StudyId = inToStudyId where StudyId = inFromStudyId;
  Delete ignore from Markers where StudyId = inFromStudyId;
  Update ignore MapMarkers set StudyId = inToStudyId where StudyId = inFromStudyId;
  Delete ignore from MapMarkers where StudyId = inFromStudyId;
  Update Pedigrees set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update ignore Positions set StudyId = inToStudyId where StudyId = inFromStudyId;
  Delete ignore from Positions where StudyId = inFromStudyId;
  Update PedigreePositions set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update ignore Regions set StudyId = inToStudyId where StudyId = inFromStudyId;
  Delete ignore from Regions where StudyId = inFromStudyId;
  Update Servers set StudyId = inToStudyId where StudyId = inFromStudyId;
  -- Models and TP2MP both reference PedPosIds, so they're already taken care of.

END;
//
Delimiter ;

DROP PROCEDURE IF EXISTS UnloadModelParts;
DELIMITER //

CREATE PROCEDURE UnloadModelParts ()
BEGIN

  DECLARE version char(64) DEFAULT '$Id: AltL_server.sql 118 2010-10-07 15:24:05Z whv001 $';

  Select Concat(
    'Insert into DModelParts (MPId, DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen) values (',
	convert(MPId,char),',',
	format(DGF,32),',',
	format(BigPen,32),',',
	format(BigLittlePen,32),',',
	format(LittleBigPen,32),',',
	format(LittlePen,32),');')
    from DModelParts;

  Select Concat(
    'Insert into QModelParts (MPId, DGF, BigPenMean, BigLittlePenMean, LittleBigPenMean, LittlePenMean, BigPenSD, BigLittlePenSD, LittleBigPenSD, LittlePenSD, Threshold) values (',
	convert(MPId,char),',',
	format(DGF,32),',',
	format(BigPenMean,32),',',
	format(BigLittlePenMean,32),',',
	format(LittleBigPenMean,32),',',
	format(LittlePenMean,32),',',
	format(BigPenSD,32),',',
	format(BigLittlePenSD,32),',',
	format(LittleBigPenSD,32),',',
	format(LittlePenSD,32),',',
	format(Threshold,32),');')
    from QModelParts;

END;
//
Delimiter ;

DROP PROCEDURE IF EXISTS Q;
DELIMITER //

CREATE PROCEDURE Q (IN inWhich varchar(16))
BEGIN

  DECLARE localTotal int;
  DECLARE localIgnore int;

WholeThing: LOOP

  -- Progress: Snapshot of statistics for the entire database
  IF inWhich = 'Progress' THEN
    Select Now(), count(*), sum(RuntimeCostSec) 'Sum', max(RuntimeCostSec) 'Max', avg(RuntimeCostSec) 'Avg'
	from Models where Likelihood is not null;
    Leave WholeThing;
  END IF;

  -- Servers: Servers-table-only status for all servers that are still running or were active in the last 4 hours
  IF inWhich = 'Servers' THEN
    Select ServerId, HostName, StudyId, PedigreeRegEx, PedigreeNotRegex, StartTime, CurrentPedPosId, CurrentLimit, LastHeartBeat, ExitStatus
	from Servers where LastHeartbeat > SubTime(Now(), '04:00') order by StudyId, ExitStatus, ServerId;
    Leave WholeThing;
  END IF;

  -- Recent: Work rate and status for servers that are still running or were active in the last hour
  IF inWhich = 'Recent' THEN
    Select s.StudyId, s.PedigreeRegEx, s.ServerId, count(*), sum(m.RuntimeCostSec) 'Sum', max(m.RuntimeCostSec) 'Max', avg(m.RuntimeCostSec) 'Avg', s.LastHeartbeat, s.ExitStatus
	from Servers s, Models m where m.ServerId = s.ServerId AND m.Likelihood is NOT NULL AND (s.LastHeartbeat > SubTime(Now(), '01:00') OR ExitStatus IS NULL)
	group by s.ServerId order by s.StudyId, s.ExitStatus, s.PedigreeRegEx;
    Leave WholeThing;
  END IF;

  -- Active: Work rate and status for servers that are still running
  IF inWhich = 'Active' THEN
    Select s.StudyId, s.PedigreeRegEx, s.ServerId, count(*), sum(m.RuntimeCostSec) 'Sum', max(m.RuntimeCostSec) 'Max', avg(m.RuntimeCostSec) 'Avg', s.LastHeartbeat
	from Servers s, Models m where m.ServerId = s.ServerId AND m.Likelihood is NOT NULL AND s.ExitStatus IS NULL
	group by s.ServerId order by s.StudyId, s.PedigreeRegEx;
    Leave WholeThing;
  END IF;

  -- Delta: Loop showing overall work progress for the last minute using servers active in the last hour
  IF inWhich = 'Delta' THEN
    MPS: LOOP
      Select SubTime(Now(), '01:00') into @startTime from dual;
      Select count(*) into localTotal from Servers s, Models m where m.ServerId = s.ServerId AND m.Likelihood is NOT NULL AND s.LastHeartbeat > @startTime;
      Select Sleep (60) into localIgnore;
      Select Now(), (count(*) - localTotal)/60.0 'MPS' from Servers s, Models m where m.ServerId = s.ServerId AND m.Likelihood is NOT NULL AND s.LastHeartbeat > @startTime;
    END LOOP MPS;
    Leave WholeThing;
  END IF;

  -- Free: Show unallocated work by StudyId/PedPosId with SingleModelRuntime
  IF inWhich = 'Free' THEN
    Select PP.StudyId, PP.PedPosId, PP.SingleModelRuntime, count(*) from
      PedigreePositions PP, Models M where
      PP.PedPosId = M.PedPosId AND M.ServerId is NULL group by PP.StudyId, PP.PedPosId;
    Leave WholeThing;
  END IF;

  IF inWhich = 'Reconcile' THEN
    Update Servers set ExitStatus = 42 where ConnectionId NOT IN (Select ID from INFORMATION_SCHEMA.PROCESSLIST) AND ExitStatus IS NULL;
    Leave WholeThing;
  END IF;

   Create temporary table Q_help (Keyword varchar(32), Action varchar(120));
   Insert into Q_help (Keyword, Action) values
	('Progress','Snapshot of statistics for the entire database'),
	('Servers', 'Servers-table-only status for all servers that are still running or were active in the last 4 hours'),
	('Recent', 'Work rate and status for servers that are still running or were active in the last hour'),
	('Active', 'Work rate and status for servers that are still running'),
	('Delta', 'Loop showing overall work progress for the last minute using servers active in the last hour'),
	('Free', 'Show unallocated work by StudyId/PedPosId with SingleModelRuntime'),
	('Reconcile', 'Mark any servers not really in processlist with ExitStatus 42');
   Select * from Q_help;
   Drop table Q_help;

  Leave WholeThing;
END LOOP WholeThing;

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS SMRT_stats;
DELIMITER //

CREATE PROCEDURE SMRT_stats (IN inLowStudyId int, IN inHighStudyId int)
BEGIN

  DECLARE localNDim int;
  DECLARE localModels int;

WholeThing: LOOP

Select a.StudyId, a.LiabilityClassCnt 'Liability Classes', a.ImprintingFlag 'Imprinting?', count(*) 'Ped/Pos count', 
  sum(b.SingleModelRuntime) 'Total SMRTs (s)', max(b.SingleModelRuntime) 'Max SMRT (s)', (4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt 'Dimensions',
  (1+
    8*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)+
    2*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)+
    4*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)+
    4*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-2))/3+
    pow(2,(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)))
  ) 'DCUHRE models',
  concat(
    convert(
      convert(
        (1+
          8*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)+
          2*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)+
          4*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)+
          4*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-2))/3+
          pow(2,(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)))
        )
        * 2 * sum(b.SingleModelRuntime) / 86400,
        unsigned
      ),
    char),
    ' ',
    convert(
      sec_to_time(
        (1+
          8*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)+
          2*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)+
          4*((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)+
          4*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-1)*(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)-2))/3+
          pow(2,(((4+(a.ImprintingFlag='y'))*a.LiabilityClassCnt)))
        ) * 2 * sum(b.SingleModelRuntime) % 86400
      ),
      char
    )
  )
  'Serial time (Days H:M:S)'
  from Studies a, SingleModelRuntimes b where a.StudyId = b.StudyId AND a.StudyId between inLowStudyId and inHighStudyId group by a.StudyId;

  Select StudyId, PedigreeSId, count(*) 'Pos count', sum(SingleModelRuntime) 'Total SMRT (s)', max(SingleModelRuntime) 'Max SMRT (s)'
	from SingleModelRuntimes where StudyId between inLowStudyId and inHighStudyId group by PedigreeSId, StudyId order by max(SingleModelRuntime) desc limit 10;

  Select 'Max SMRT values less than 10 seconds can be Polynomial runs (in 16G), but if any are over 20 seconds they must be NonPolynomial runs' as 'Note:';

  Leave WholeThing;
END LOOP WholeThing;

END
//
DELIMITER ;

