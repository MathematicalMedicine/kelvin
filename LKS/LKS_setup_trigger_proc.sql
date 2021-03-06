-- Copyright (C) 2022 Mathematical Medicine LLC
-- 
-- This program is free software: you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the Free
-- Software Foundation, either version 3 of the License, or (at your option)
-- any later version.
-- 
-- You should have received a copy of the GNU General Public License along
-- with this program. If not, see <https://www.gnu.org/licenses/>.

--
-- Dangerous @-sign GLOBALS...look for test-set-act sequences that can deadlock and break during the act
-- because even if they are in a transaction, the set is not something that gets rolled back.
--
-- @outMarkerSetId is...
--   1. Set to LAST_INSERT_ID into MarkerSetLikelihood at the very end of GetMarkerSetLikelihood
--   2. Set to outLC2MPId or -1 at 3 places in GetWork, eeeeeuuuugh!
--   3. Set to LAST_INSERT_ID into MarkerSetLikelihood in CacheMarkerWork
--   4. In a diagnostic "print" and as key in select/counting, inserting and updating MarkerSetLikelihood_MCMC in PutWork
--
-- @MCMC_flag is...
--   1. Set to 1 or 0 in GetWork
--   2. Tested in CacheCombinedWork
--   3. Tested in PutWork
--
-- @outLCM1MPId is...
--   1. Conditionally set to -1 in GetWork if realLocusListType = 1 (right after call to CacheMarkerWork)
--   2. Tested for >= 0 in PutWork
--   AND NEVER CHANGES?! This has _got_ to be broken.
--
-- I wonder about @outModelId...but it is local to PutWork
--

-- NOTE: A number of comments are marked as follows:
    -- DEBUG_DIAGTABLE 
-- These are for statements that populate a "Diag" table that logs a whole
-- lot of diagnostic output. It gets really really large really really quickly,
-- tho, so we generally will not have it turned on.
-- 
-- Others are marked as follows:
    -- DEBUG_RMTABLE 
-- This is for the "RegionModels" table, which is a much smaller diagnostic
-- tool.
--
-- These are being turned "on" and "off" by search-and-replacing those strings
-- (without leading whitespace) and either appending or removing a trailing
-- newline.


-- All the careful tracking of available work is no longer necessary given our pipelined approach. It has been
-- being phased out for some time, and finally we're dropping it and instead going to use FreeModels as a flag
-- that is set in the pipeline after client and server-set runs.

DROP TRIGGER IF EXISTS BeforeInsertDiag;

-- DEBUG_DIAGTABLE DELIMITER //
-- DEBUG_DIAGTABLE CREATE TRIGGER BeforeInsertDiag BEFORE INSERT ON Diag FOR EACH ROW
-- DEBUG_DIAGTABLE BEGIN
-- DEBUG_DIAGTABLE   IF new.ConnectionId IS NULL
-- DEBUG_DIAGTABLE   THEN
-- DEBUG_DIAGTABLE     SET new.ConnectionId = CONNECTION_ID();
-- DEBUG_DIAGTABLE   END IF;
-- DEBUG_DIAGTABLE END;
-- DEBUG_DIAGTABLE //
-- DEBUG_DIAGTABLE DELIMITER ;

-- DELIMITER ;
DROP TRIGGER IF EXISTS DecMarkers;

-- DELIMITER //
-- CREATE TRIGGER DecMarkers AFTER UPDATE ON MarkerSetLikelihood FOR EACH ROW
-- BEGIN
--   Update PedigreePositions set FreeModels = FreeModels - 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.ServerId IS NULL AND
--     NEW.ServerId IS NOT NULL;
--   Update PedigreePositions set FreeModels = FreeModels + 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.ServerId IS NOT NULL AND
--     NEW.ServerId IS NULL;
--   Update PedigreePositions set PendingLikelihoods = PendingLikelihoods - 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.Likelihood IS NULL AND
--     NEW.Likelihood IS NOT NULL;
--   Update PedigreePositions set PendingLikelihoods = PendingLikelihoods + 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.Likelihood IS NOT NULL AND
--     NEW.Likelihood IS NULL;
-- END;
-- //
-- DELIMITER ;

DROP TRIGGER IF EXISTS IncMarkers;

-- DELIMITER //
-- CREATE TRIGGER IncMarkers AFTER INSERT ON MarkerSetLikelihood FOR EACH ROW
-- BEGIN
--   Update PedigreePositions set FreeModels = FreeModels + 1 where
--     PedPosId = NEW.PedPosId AND
--     NEW.ServerId IS NULL AND
--     NEW.Likelihood IS NULL;
--   Update PedigreePositions set PendingLikelihoods = PendingLikelihoods + 1 where
--     PedPosId = NEW.PedPosId AND
--     NEW.Likelihood IS NULL;
-- END;
-- //
-- DELIMITER ;

DROP TRIGGER IF EXISTS DecModels;

-- DELIMITER //
-- CREATE TRIGGER DecModels AFTER UPDATE ON Models FOR EACH ROW
-- BEGIN
--   Update PedigreePositions set FreeModels = FreeModels - 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.ServerId IS NULL AND
--     NEW.ServerId IS NOT NULL;
--   Update PedigreePositions set FreeModels = FreeModels + 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.ServerId IS NOT NULL AND
--     NEW.ServerId IS NULL;
--   Update PedigreePositions set PendingLikelihoods = PendingLikelihoods - 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.Likelihood IS NULL AND
--     NEW.Likelihood IS NOT NULL;
--   Update PedigreePositions set PendingLikelihoods = PendingLikelihoods + 1 where
--     PedPosId = NEW.PedPosId AND
--     OLD.Likelihood IS NOT NULL AND
--     NEW.Likelihood IS NULL;
-- END;
-- //
-- DELIMITER ;

DROP TRIGGER IF EXISTS IncModels;

-- DELIMITER //
-- CREATE TRIGGER IncModels AFTER INSERT ON Models FOR EACH ROW
-- BEGIN
--   Update PedigreePositions set FreeModels = FreeModels + 1 where
--     PedPosId = NEW.PedPosId AND
--     NEW.ServerId IS NULL AND
--     NEW.Likelihood IS NULL;
--   Update PedigreePositions set PendingLikelihoods = PendingLikelihoods + 1 where
--     PedPosId = NEW.PedPosId AND
--     NEW.Likelihood IS NULL;
-- END;
-- //
-- DELIMITER ;

DROP PROCEDURE IF EXISTS BadScaling;

DELIMITER //
-- BadScaling is called after the addition of new trait positions or new pedigree genotype
-- maps. It (re)populates the MapMarkers scaling information and PedigreePositions PedTraitPosCM
-- column. It uses a bad linear interpolation approach for marker-bounded positions, and a
-- smarter extrapolation approach for the rest.

CREATE PROCEDURE BadScaling (IN inStudyId int)
BEGIN
  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE no_rows_indicator INT DEFAULT 0;
  DECLARE localMapId INT;
  DECLARE localMarkerName char(16);
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('BadScaling: called w/ ', convert(IFNULL(inStudyId,'NULL'),char)));

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
    Select MapId, MarkerName into localMapId, localMarkerName from MapMarkers where
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

  Delete from MapMarkers where MarkerName = 'Dummy' AND StudyId = inStudyId;

  -- Yeah, this is just a diagnostic, a way of bailing early to test incrementally.
  Leave WholeThing;

END LOOP WholeThing;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('BadScaling: returning.'));

END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetMarkerSetLikelihood;
DELIMITER //
CREATE PROCEDURE GetMarkerSetLikelihood (
 IN inPedPosId int, 
 IN inAnalysisId int,
 IN inRegionNo int, IN inParentRegionNo int, IN inParentRegionError real, IN inParentRegionSplitDir int,
 OUT outRegionId int, OUT outMarkerCount int, OUT outLikelihood real
)
BEGIN
    DECLARE version char(96) DEFAULT '$Id$';
    DECLARE WNE_indicator INT DEFAULT 0; -- For when we know the Likelihood _will_ not exist
    DECLARE outMarkerSetId INT DEFAULT -1;
    DECLARE no_rows_indicator INT DEFAULT 0;
    DECLARE localStudyId INT;
    DECLARE localChromosomeNo INT;
    DECLARE localRefTraitPosCM DOUBLE;
    DECLARE no_rows CONDITION FOR 1329;
    DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetMarkerSetLikelihood: called w/ ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inPedPosId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inAnalysisId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inRegionNo,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionNo,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionError,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionSplitDir,'NULL'),char)
    -- DEBUG_DIAGTABLE ));

    Start transaction;

    -- Get the RegionId 

    Select StudyId, ChromosomeNo, RefTraitPosCM INTO localStudyId, localChromosomeNo, localRefTraitPosCM
    from PedigreePositions where PedPosId = inPedPosId;

    Set no_rows_indicator = 0;
    Select RegionId into outRegionId from Regions where
    StudyId = localStudyId AND AnalysisId = inAnalysisId AND ChromosomeNo = localChromosomeNo AND 
    RefTraitPosCM = localRefTraitPosCM AND RegionNo = inRegionNo;

    If no_rows_indicator THEN
      Insert into Regions (StudyId, AnalysisId, ChromosomeNo, RefTraitPosCM, RegionNo, ParentRegionNo, ParentRegionError, ParentRegionSplitDir) values
      (localStudyId, inAnalysisId, localChromosomeNo, localRefTraitPosCM, inRegionNo, inParentRegionNo, inParentRegionError, inParentRegionSplitDir);
      Select LAST_INSERT_ID() INTO outRegionId;
    END IF;

    SET no_rows_indicator = 0;
    Select MarkerSetId, Likelihood, MarkerCount into outMarkerSetId, outLikelihood, outMarkerCount from MarkerSetLikelihood 
	where PedPosId = inPedPosId order by MarkerCount desc limit 1;
    -- insert a new entry in MarkerSetLikelihood
    IF no_rows_indicator THEN
      Insert into MarkerSetLikelihood (PedPosId) values  (inPedPosId);
      select LAST_INSERT_ID() into outMarkerSetId;
    END IF;
    set @outMarkerSetId = outMarkerSetId;
    commit;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetMarkerSetLikelihood: returning w/ ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(outRegionId,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(outMarkerCount,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(outLikelihood,'NULL'),char),', '
  -- DEBUG_DIAGTABLE   ));

END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetAnalysisId;
DELIMITER //
CREATE PROCEDURE GetAnalysisId (
 IN inStudyId int, 
 IN inPedigreeRegEx varchar(1025),
 IN inPedigreeNotRegEx varchar(1025),
 OUT outAnalysisId int)
BEGIN
    DECLARE version char(96) DEFAULT '$Id$';
    -- This now relies exclusively on the unique key, which really should be PedigreRegEx and PedigreeNotRegEx, but
    -- those have become unreasonably long.

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetAnalysisId: called w/ ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inStudyId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inPedigreeRegEx,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inPedigreeNotRegEx,'NULL'),char)));

    Insert ignore into Analyses (StudyId, Uniquey, PedigreeRegEx, PedigreeNotRegEx) values
      (inStudyId, COMPRESS(CONCAT(inPedigreeRegEx,'/',inPedigreeNotRegEx)), inPedigreeRegEx, inPedigreeNotRegEx);

    Select AnalysisId INTO outAnalysisId
    from Analyses where StudyId = inStudyId AND PedigreeRegEx = inPedigreeRegEx AND PedigreeNotRegEx = inPedigreeNotRegEx;

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetAnalysisId: returning w/ ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(outAnalysisId,'NULL'),char)));
      
END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetDLikelihood;
DELIMITER //

-- GetLikelihood is called by the kelvin client to request and retrieve Likelihoods. All possible fixed-model
-- parameters are passed in the call. They should be -1 if they're N/A. If a corresponding Likelihood
-- is found, it is returned, otherwise the call is treated as an asynchronous request and a NULL
-- is returned for the Likelihood.

CREATE PROCEDURE GetDLikelihood (
 IN inPedPosId int, IN inDGF real, 
 IN inLC1BigPen real, IN inLC1BigLittlePen real, IN inLC1LittleBigPen real, IN inLC1LittlePen real,
 IN inLC2BigPen real, IN inLC2BigLittlePen real, IN inLC2LittleBigPen real, IN inLC2LittlePen real,
 IN inLC3BigPen real, IN inLC3BigLittlePen real, IN inLC3LittleBigPen real, IN inLC3LittlePen real,
 IN inAnalysisId int,
 IN inRegionNo int, IN inParentRegionNo int, IN inParentRegionError real, IN inParentRegionSplitDir int,
 OUT outRegionId int, OUT outMarkerCount int, OUT outLikelihood real
)
BEGIN
    DECLARE version char(96) DEFAULT '$Id$';
    DECLARE WNE_indicator INT DEFAULT 0; -- For when we know the Likelihood will not exist
    DECLARE localStudyId INT;
    DECLARE localChromosomeNo INT;
    DECLARE localRefTraitPosCM DOUBLE;
    DECLARE localLC1MPId INT;
    DECLARE localLC2MPId INT;
    DECLARE localLC3MPId INT;
    DECLARE localModelId INT;
    DECLARE no_rows_indicator INT DEFAULT 0;
    DECLARE no_rows CONDITION FOR 1329;
    DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetDLikelihood: called w/ ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inPedPosId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inDGF,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC1BigPen,'NULL'),char),', ', convert(IFNULL(inLC1BigLittlePen,'NULL'),char),', ', convert(IFNULL(inLC1LittleBigPen,'NULL'),char),', ', convert(IFNULL(inLC1LittlePen,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC2BigPen,'NULL'),char),', ', convert(IFNULL(inLC2BigLittlePen,'NULL'),char),', ', convert(IFNULL(inLC2LittleBigPen,'NULL'),char),', ', convert(IFNULL(inLC2LittlePen,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC3BigPen,'NULL'),char),', ', convert(IFNULL(inLC3BigLittlePen,'NULL'),char),', ', convert(IFNULL(inLC3LittleBigPen,'NULL'),char),', ', convert(IFNULL(inLC3LittlePen,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inAnalysisId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inRegionNo,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionNo,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionError,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionSplitDir,'NULL'),char)
    -- DEBUG_DIAGTABLE   ));

    Start transaction;

    -- Get the RegionId 

    Select StudyId, ChromosomeNo, RefTraitPosCM INTO localStudyId, localChromosomeNo, localRefTraitPosCM
    from PedigreePositions where PedPosId = inPedPosId;

    Set no_rows_indicator = 0;
    Select RegionId into outRegionId from Regions where
    StudyId = localStudyId AND AnalysisId = inAnalysisId AND ChromosomeNo = localChromosomeNo AND
    RefTraitPosCM = localRefTraitPosCM AND RegionNo = inRegionNo;

    If no_rows_indicator THEN
      Insert into Regions (StudyId, AnalysisId, ChromosomeNo, RefTraitPosCM, RegionNo, ParentRegionNo, ParentRegionError, ParentRegionSplitDir) values
      (localStudyId, inAnalysisId, localChromosomeNo, localRefTraitPosCM, inRegionNo, inParentRegionNo, inParentRegionError, inParentRegionSplitDir);
      Select LAST_INSERT_ID() INTO outRegionId;
    END IF;

    -- Handle marker set likelihood (which doesn't use the trait model) separately
    IF inDGF = -1 THEN
      -- This is a marker set likelihood
      SET no_rows_indicator = 0;
      Select Likelihood, MarkerCount into outLikelihood, outMarkerCount from
	MarkerSetLikelihood where PedPosId = inPedPosId order by MarkerCount desc limit 1;
      -- insert a new entry in MarkerSetLikelihood
      IF no_rows_indicator THEN
        Insert into MarkerSetLikelihood (PedPosId) values
	  (inPedPosId);
      END IF;
    ELSE
      -- This is a trait or combined likelihood, so the trait model is relevant
      -- There's always at least one liability class...
      SET no_rows_indicator = 0;
      -- PERFORMANCE-WISE, this combination of constraint columns is explicitly indexed by a UNIQUE KEY
      Select MPId into localLC1MPId from DModelParts where
      DGF = convert(IFNULL(inDGF,'NULL'), DECIMAL(32,30)) AND BigPen = convert(IFNULL(inLC1BigPen,'NULL'), DECIMAL(32,30)) AND BigLittlePen = convert(IFNULL(inLC1BigLittlePen,'NULL'), DECIMAL(32,30)) AND 
      LittleBigPen = convert(IFNULL(inLC1LittleBigPen,'NULL'), DECIMAL(32,30)) AND LittlePen = convert(IFNULL(inLC1LittlePen,'NULL'), DECIMAL(32,30));

      IF no_rows_indicator THEN
        Insert into DModelParts (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)
        values (inDGF, inLC1BigPen, inLC1BigLittlePen, inLC1LittleBigPen, inLC1LittlePen);
        Select LAST_INSERT_ID() INTO localLC1MPId;
        SET WNE_indicator = 1;
      END IF;

      SET no_rows_indicator = 0;
      -- PERFORMANCE-WISE, this is explicitly indexed
      Select MPId into localLC2MPId from DModelParts where
      DGF = convert(IFNULL(inDGF,'NULL'), DECIMAL(32,30)) AND BigPen = convert(IFNULL(inLC2BigPen,'NULL'), DECIMAL(32,30)) AND BigLittlePen = convert(IFNULL(inLC2BigLittlePen,'NULL'), DECIMAL(32,30)) AND 
      LittleBigPen = convert(IFNULL(inLC2LittleBigPen,'NULL'), DECIMAL(32,30)) AND LittlePen = convert(IFNULL(inLC2LittlePen,'NULL'), DECIMAL(32,30));
    
      IF no_rows_indicator THEN
        Insert into DModelParts (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)
        values (inDGF, inLC2BigPen, inLC2BigLittlePen, inLC2LittleBigPen, inLC2LittlePen);
        Select LAST_INSERT_ID() INTO localLC2MPId;
        SET WNE_indicator = 1;
      END IF;

      SET no_rows_indicator = 0;
      -- PERFORMANCE-WISE, this is explicitly indexed
      Select MPId into localLC3MPId from DModelParts where
      DGF = convert(IFNULL(inDGF,'NULL'), DECIMAL(32,30)) AND BigPen = convert(IFNULL(inLC3BigPen,'NULL'), DECIMAL(32,30)) AND BigLittlePen = convert(IFNULL(inLC3BigLittlePen,'NULL'), DECIMAL(32,30)) AND 
      LittleBigPen = convert(IFNULL(inLC3LittleBigPen,'NULL'), DECIMAL(32,30)) AND LittlePen = convert(IFNULL(inLC3LittlePen,'NULL'), DECIMAL(32,30));
    
      IF no_rows_indicator THEN
        Insert into DModelParts (DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen)
        values (inDGF, inLC3BigPen, inLC3BigLittlePen, inLC3LittleBigPen, inLC3LittlePen);
        Select LAST_INSERT_ID() INTO localLC3MPId;
        SET WNE_indicator = 1;
      END IF;

      -- Maybe attempt to insert a row in Models table
      IF WNE_indicator THEN
        -- Some of the model parts are not even in the table, so this model was not either,
	-- although that does not guarentee that it isn't there now.
        Insert ignore into Models (PedPosId, LC1MPId, LC2MPId, LC3MPId) values
  	(inPedPosId, localLC1MPId, localLC2MPId, localLC3MPId);
      END IF;

      -- Go for the Likelihood! If it doesn't exist, request it and return the NULL.
      SET no_rows_indicator = 0;
      -- PERFORMANCE-WISE, this is a PK query, therefore indexed
      Select Likelihood, MarkerCount into outLikelihood, outMarkerCount from
      Models where PedPosId = inPedPosId AND LC1MPId = localLC1MPId AND
        LC2MPId = localLC2MPId AND LC3MPId = localLC3MPId order by MarkerCount desc limit 1;

      IF no_rows_indicator THEN
        -- insert ignore because we just might have inserted it earlier if the
        -- WNE_indicator was set
        Insert ignore into Models (PedPosId, LC1MPId, LC2MPId, LC3MPId) values
        (inPedPosId, localLC1MPId, localLC2MPId, localLC3MPId);
        Select LAST_INSERT_ID() INTO localModelId;
        -- DEBUG_RMTABLE Insert ignore into RegionModels (RegionId, ModelId) values (outRegionId, localModelId); -- still a valuable diagnostic
      END IF;
    END IF;
    Commit;

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetDLikelihood: returning w/ ',
    -- DEBUG_DIAGTABLE convert(IFNULL(outRegionId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE convert(IFNULL(outMarkerCount,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE convert(IFNULL(outLikelihood,'NULL'),char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetQLikelihood;
DELIMITER //

-- GetLikelihood is called by the kelvin client to request and retrieve Likelihoods. All possible fixed-model
-- parameters are passed in the call. They should be -1 if they're N/A. If a corresponding Likelihood
-- is found, it is returned, otherwise the call is treated as an asynchronous request and a NULL
-- is returned for the Likelihood.

CREATE PROCEDURE GetQLikelihood (
 IN inPedPosId int, IN inDGF real, 
 IN inLC1BigMean real, IN inLC1BigLittleMean real, IN inLC1LittleBigMean real, IN inLC1LittleMean real,
 IN inLC2BigMean real, IN inLC2BigLittleMean real, IN inLC2LittleBigMean real, IN inLC2LittleMean real,
 IN inLC3BigMean real, IN inLC3BigLittleMean real, IN inLC3LittleBigMean real, IN inLC3LittleMean real,
 IN inLC1BigSD real, IN inLC1BigLittleSD real, IN inLC1LittleBigSD real, IN inLC1LittleSD real,
 IN inLC2BigSD real, IN inLC2BigLittleSD real, IN inLC2LittleBigSD real, IN inLC2LittleSD real,
 IN inLC3BigSD real, IN inLC3BigLittleSD real, IN inLC3LittleBigSD real, IN inLC3LittleSD real,
 IN inLC1Threshold real, IN inLC2Threshold real, IN inLC3Threshold real,
 IN inAnalysisId int,
 IN inRegionNo int, IN inParentRegionNo int, IN inParentRegionError real, IN inParentRegionSplitDir int,
 OUT outRegionId int, OUT outMarkerCount int, OUT outLikelihood real
)
BEGIN
    DECLARE version char(96) DEFAULT '$Id$';
    DECLARE WNE_indicator INT DEFAULT 0; -- For when we know the Likelihood _will_ not exist
    DECLARE localStudyId INT;
    DECLARE localChromosomeNo INT;
    DECLARE localRefTraitPosCM DOUBLE;
    DECLARE localLC1MPId INT;
    DECLARE localLC2MPId INT;
    DECLARE localLC3MPId INT;
    DECLARE localModelId INT;
    DECLARE no_rows_indicator INT DEFAULT 0;
    DECLARE no_rows CONDITION FOR 1329;
    DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetQLikelihood: called w/ ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inPedPosId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inDGF,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC1BigMean,'NULL'),char),', ', convert(IFNULL(inLC1BigLittleMean,'NULL'),char),', ', convert(IFNULL(inLC1LittleBigMean,'NULL'),char),', ', convert(IFNULL(inLC1LittleMean,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC2BigMean,'NULL'),char),', ', convert(IFNULL(inLC2BigLittleMean,'NULL'),char),', ', convert(IFNULL(inLC2LittleBigMean,'NULL'),char),', ', convert(IFNULL(inLC2LittleMean,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC3BigMean,'NULL'),char),', ', convert(IFNULL(inLC3BigLittleMean,'NULL'),char),', ', convert(IFNULL(inLC3LittleBigMean,'NULL'),char),', ', convert(IFNULL(inLC3LittleMean,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC1BigSD,'NULL'),char),', ', convert(IFNULL(inLC1BigLittleSD,'NULL'),char),', ', convert(IFNULL(inLC1LittleBigSD,'NULL'),char),', ', convert(IFNULL(inLC1LittleSD,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC2BigSD,'NULL'),char),', ', convert(IFNULL(inLC2BigLittleSD,'NULL'),char),', ', convert(IFNULL(inLC2LittleBigSD,'NULL'),char),', ', convert(IFNULL(inLC2LittleSD,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC3BigSD,'NULL'),char),', ', convert(IFNULL(inLC3BigLittleSD,'NULL'),char),', ', convert(IFNULL(inLC3LittleBigSD,'NULL'),char),', ', convert(IFNULL(inLC3LittleSD,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inLC1Threshold,'NULL'),char),', ', convert(IFNULL(inLC2Threshold,'NULL'),char),', ', convert(IFNULL(inLC3Threshold,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inAnalysisId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inRegionNo,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionNo,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionError,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE   convert(IFNULL(inParentRegionSplitDir,'NULL'),char)
    -- DEBUG_DIAGTABLE   ));

    Start transaction;

    -- Get the RegionId

    Select StudyId, ChromosomeNo, RefTraitPosCM INTO localStudyId, localChromosomeNo, localRefTraitPosCM
    from PedigreePositions where PedPosId = inPedPosId;

    Set no_rows_indicator = 0;
    Select RegionId into outRegionId from Regions where
    StudyId = localStudyId AND AnalysisId = inAnalysisId AND ChromosomeNo = localChromosomeNo AND 
    RefTraitPosCM = localRefTraitPosCM AND RegionNo = inRegionNo;

    If no_rows_indicator THEN
      Insert into Regions (StudyId, AnalysisId, ChromosomeNo, RefTraitPosCM, RegionNo, ParentRegionNo, ParentRegionError, ParentRegionSplitDir) values
      (localStudyId, inAnalysisId, localChromosomeNo, localRefTraitPosCM, inRegionNo, inParentRegionNo, inParentRegionError, inParentRegionSplitDir);
      Select LAST_INSERT_ID() INTO outRegionId;
    END IF;

    -- Handle marker set likelihood (which doesn't use the trait model) separately
    IF inDGF = -1 THEN
      -- This is a marker set likelihood
      SET no_rows_indicator = 0;
      Select Likelihood, MarkerCount into outLikelihood, outMarkerCount from
	MarkerSetLikelihood where PedPosId = inPedPosId order by MarkerCount desc limit 1;
      -- insert a new entry in MarkerSetLikelihood
      IF no_rows_indicator THEN
        Insert into MarkerSetLikelihood (PedPosId) values
	  (inPedPosId);
      END IF;
    ELSE
      -- This is a trait or combined likelihood, so the trait model is relevant
      -- There's always at least one liability class...
      SET no_rows_indicator = 0;
      Select MPId into localLC1MPId from QModelParts where
        DGF = convert(IFNULL(inDGF,'NULL'), DECIMAL(32,30)) AND BigMean = convert(IFNULL(inLC1BigMean,'NULL'), DECIMAL(32,30)) AND BigLittleMean = convert(IFNULL(inLC1BigLittleMean,'NULL'), DECIMAL(32,30)) AND 
        LittleBigMean = convert(IFNULL(inLC1LittleBigMean,'NULL'), DECIMAL(32,30)) AND LittleMean = convert(IFNULL(inLC1LittleMean,'NULL'), DECIMAL(32,30)) AND BigSD=convert(IFNULL(inLC1BigSD,'NULL'), Decimal(32,30)) AND BigLittleSD=convert(IFNULL(inLC1BigLittleSD,'NULL'), Decimal(32,30)) AND LittleBigSD=convert(IFNULL(inLC1LittleBigSD,'NULL'), Decimal(32,30)) AND LittleSD=convert(IFNULL(inLC1LittleSD,'NULL'), Decimal(32,30)) AND Threshold=convert(IFNULL(inLC1Threshold,'NULL'), Decimal(32,30));

      IF no_rows_indicator THEN
        Insert into QModelParts (DGF, BigMean, BigLittleMean, LittleBigMean, LittleMean, BigSD, BigLittleSD, LittleBigSD, LittleSD, Threshold)
        values (inDGF, inLC1BigMean, inLC1BigLittleMean, inLC1LittleBigMean, inLC1LittleMean, inLC1BigSD, inLC1BigLittleSD, inLC1LittleBigSD, inLC1LittleSD, inLC1Threshold);
        Select LAST_INSERT_ID() INTO localLC1MPId;
        SET WNE_indicator = 1;
      END IF;

      SET no_rows_indicator = 0;
      Select MPId into localLC2MPId from QModelParts where
      DGF = convert(IFNULL(inDGF,'NULL'), DECIMAL(32,30)) AND BigMean = convert(IFNULL(inLC2BigMean,'NULL'), DECIMAL(32,30)) AND BigLittleMean = convert(IFNULL(inLC2BigLittleMean,'NULL'), DECIMAL(32,30)) AND 
      LittleBigMean = convert(IFNULL(inLC2LittleBigMean,'NULL'), DECIMAL(32,30)) AND LittleMean = convert(IFNULL(inLC2LittleMean,'NULL'), DECIMAL(32,30)) AND BigSD=convert(IFNULL(inLC2BigSD,'NULL'), Decimal(32,30)) AND BigLittleSD=convert(IFNULL(inLC2BigLittleSD,'NULL'), Decimal(32,30)) AND LittleBigSD=convert(IFNULL(inLC2LittleBigSD,'NULL'), Decimal(32,30)) AND LittleSD=convert(IFNULL(inLC2LittleSD,'NULL'), Decimal(32,30)) AND Threshold=convert(IFNULL(inLC2Threshold,'NULL'), Decimal(32,30));
    
      IF no_rows_indicator THEN
        Insert into QModelParts (DGF, BigMean, BigLittleMean, LittleBigMean, LittleMean, BigSD, BigLittleSD, LittleBigSD, LittleSD, Threshold)
        values (inDGF, inLC2BigMean, inLC2BigLittleMean, inLC2LittleBigMean, inLC2LittleMean, inLC2BigSD, inLC2BigLittleSD, inLC2LittleBigSD, inLC2LittleSD, inLC2Threshold);
        Select LAST_INSERT_ID() INTO localLC2MPId;
        SET WNE_indicator = 1;
      END IF;

      SET no_rows_indicator = 0;
      Select MPId into localLC3MPId from QModelParts where
      DGF = convert(IFNULL(inDGF,'NULL'), DECIMAL(32,30)) AND BigMean = convert(IFNULL(inLC3BigMean,'NULL'), DECIMAL(32,30)) AND BigLittleMean = convert(IFNULL(inLC3BigLittleMean,'NULL'), DECIMAL(32,30)) AND 
      LittleBigMean = convert(IFNULL(inLC3LittleBigMean,'NULL'), DECIMAL(32,30)) AND LittleMean = convert(IFNULL(inLC3LittleMean,'NULL'), DECIMAL(32,30)) AND BigSD=convert(IFNULL(inLC3BigSD,'NULL'), Decimal(32,30)) AND BigLittleSD=convert(IFNULL(inLC3BigLittleSD,'NULL'), Decimal(32,30)) AND LittleBigSD=convert(IFNULL(inLC3LittleBigSD,'NULL'), Decimal(32,30)) AND LittleSD=convert(IFNULL(inLC3LittleSD,'NULL'), Decimal(32,30)) AND Threshold=convert(IFNULL(inLC3Threshold,'NULL'), Decimal(32,30));
    
      IF no_rows_indicator THEN
        Insert into QModelParts (DGF, BigMean, BigLittleMean, LittleBigMean, LittleMean, BigSD, BigLittleSD, LittleBigSD, LittleSD, Threshold)
        values (inDGF, inLC3BigMean, inLC3BigLittleMean, inLC3LittleBigMean, inLC3LittleMean, inLC3BigSD, inLC3BigLittleSD, inLC3LittleBigSD, inLC3LittleSD, inLC3Threshold);
        Select LAST_INSERT_ID() INTO localLC3MPId;
        SET WNE_indicator = 1;
      END IF;

      -- Maybe attempt to insert a row in Models table
      IF WNE_indicator THEN
        -- Some of the model parts are not even in the table, so this model was not either,
	-- although that does not guarentee that it isn't there now.
        Insert ignore into Models (PedPosId, LC1MPId, LC2MPId, LC3MPId) values
  	(inPedPosId, localLC1MPId, localLC2MPId, localLC3MPId);
        Select LAST_INSERT_ID() INTO localModelId;
      END IF;

      -- Go for the Likelihood! If it doesn't exist, request it and return the NULL.
      SET no_rows_indicator = 0;
      Select Likelihood, MarkerCount into outLikelihood, outMarkerCount from
      Models where PedPosId = inPedPosId AND LC1MPId = localLC1MPId AND
        LC2MPId = localLC2MPId AND LC3MPId = localLC3MPId order by MarkerCount desc limit 1;
	
      IF no_rows_indicator THEN
        -- insert ignore because we just might have inserted it earlier if the
        -- WNE_indicator was set
        Insert ignore into Models (PedPosId, LC1MPId, LC2MPId, LC3MPId) values
        (inPedPosId, localLC1MPId, localLC2MPId, localLC3MPId);
        Select LAST_INSERT_ID() INTO localModelId;
        -- DEBUG_RMTABLE Insert ignore into RegionModels (RegionId, ModelId) values (outRegionId, localModelId); -- still a valuable diagnostic
      END IF;
    END IF;
    Commit;

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetQLikelihood: returning w/ ',
    -- DEBUG_DIAGTABLE convert(IFNULL(outRegionId,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE convert(IFNULL(outMarkerCount,'NULL'),char),', ',
    -- DEBUG_DIAGTABLE convert(IFNULL(outLikelihood,'NULL'),char)));

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
  DECLARE version char(96) DEFAULT '$Id$';

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetDParts: called w/ ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(inLC1MPId,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(inLC2MPId,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(inLC3MPId,'NULL'),char)
  -- DEBUG_DIAGTABLE   ));

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
    
  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetDParts: returning ', convert(IFNULL(outDGF,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE 	convert(IFNULL(outLC1BP,'NULL'),char),', ', convert(IFNULL(outLC1BLP,'NULL'),char),', ', convert(IFNULL(outLC1LBP,'NULL'),char),', ', convert(IFNULL(outLC1LP,'NULL'),char),', ', 
  -- DEBUG_DIAGTABLE 	convert(IFNULL(outLC2BP,'NULL'),char),', ', convert(IFNULL(outLC2BLP,'NULL'),char),', ', convert(IFNULL(outLC2LBP,'NULL'),char),', ', convert(IFNULL(outLC2LP,'NULL'),char),', ', 
  -- DEBUG_DIAGTABLE 	convert(IFNULL(outLC3BP,'NULL'),char),', ', convert(IFNULL(outLC3BLP,'NULL'),char),', ', convert(IFNULL(outLC3LBP,'NULL'),char),', ', convert(IFNULL(outLC3LP,'NULL'),char)));

END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetQParts;
DELIMITER //

CREATE PROCEDURE GetQParts (
  IN inLC1MPId int, inLC2MPId int, inLC3MPId int,
  OUT outDGF real,
  OUT outLC1BMean real, OUT outLC1BLMean real, OUT outLC1LBMean real, OUT outLC1LMean real,
  OUT outLC2BMean real, OUT outLC2BLMean real, OUT outLC2LBMean real, OUT outLC2LMean real,
  OUT outLC3BMean real, OUT outLC3BLMean real, OUT outLC3LBMean real, OUT outLC3LMean real,
  OUT outLC1BSd real, OUT outLC1BLSd real, OUT outLC1LBSd real, OUT outLC1LSd real,
  OUT outLC2BSd real, OUT outLC2BLSd real, OUT outLC2LBSd real, OUT outLC2LSd real,
  OUT outLC3BSd real, OUT outLC3BLSd real, OUT outLC3LBSd real, OUT outLC3LSd real,
  OUT outLC1Threshold real, OUT outLC2Threshold real, OUT outLC3Threshold real
)
BEGIN
  DECLARE version char(96) DEFAULT '$Id$';

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetQParts: called w/ ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(inLC1MPId,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(inLC2MPId,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE   convert(IFNULL(inLC3MPId,'NULL'),char)
  -- DEBUG_DIAGTABLE   ));

  -- Better to return the results as a result set than out parameters (which suck!)
  Select LC1.DGF,
	LC1.BigMean, LC1.BigLittleMean, LC1.LittleBigMean, LC1.LittleMean,
	LC2.BigMean, LC2.BigLittleMean, LC2.LittleBigMean, LC2.LittleMean,
	LC3.BigMean, LC3.BigLittleMean, LC3.LittleBigMean, LC3.LittleMean,
	LC1.BigSD, LC1.BigLittleSD, LC1.LittleBigSD, LC1.LittleSD,
	LC2.BigSD, LC2.BigLittleSD, LC2.LittleBigSD, LC2.LittleSD,
	LC3.BigSD, LC3.BigLittleSD, LC3.LittleBigSD, LC3.LittleSD,
        LC1.Threshold, LC2.Threshold, LC3.Threshold
  into
	outDGF,
	outLC1BMean, outLC1BLMean, outLC1LBMean, outLC1LMean, 
	outLC2BMean, outLC2BLMean, outLC2LBMean, outLC2LMean, 
	outLC3BMean, outLC3BLMean, outLC3LBMean, outLC3LMean, 
	outLC1BSD, outLC1BLSD, outLC1LBSD, outLC1LSD, 
	outLC2BSD, outLC2BLSD, outLC2LBSD, outLC2LSD, 
	outLC3BSD, outLC3BLSD, outLC3LBSD, outLC3LSD, 
        outLC1Threshold, outLC2Threshold, outLC3Threshold
  from QModelParts LC1, QModelParts LC2, QModelParts LC3
  where LC1.MPId = inLC1MPId AND LC2.MPId = inLC2MPId AND LC3.MPId = inLC3MPId;
    
  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetQParts: returning ', convert(IFNULL(outDGF,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE 	convert(IFNULL(outLC1BMean,'NULL'),char),', ', convert(IFNULL(outLC1BLMean,'NULL'),char),', ', convert(IFNULL(outLC1LBMean,'NULL'),char),', ', convert(IFNULL(outLC1LMean,'NULL'),char),', ', 
  -- DEBUG_DIAGTABLE 	convert(IFNULL(outLC2BMean,'NULL'),char),', ', convert(IFNULL(outLC2BLMean,'NULL'),char),', ', convert(IFNULL(outLC2LBMean,'NULL'),char),', ', convert(IFNULL(outLC2LMean,'NULL'),char),', ', 
  -- DEBUG_DIAGTABLE 	convert(IFNULL(outLC3BMean,'NULL'),char),', ', convert(IFNULL(outLC3BLMean,'NULL'),char),', ', convert(IFNULL(outLC3LBMean,'NULL'),char),', ', convert(IFNULL(outLC3LMean,'NULL'),char)));

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

  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE realLocusListType INT DEFAULT 0;
  DECLARE totalSampleCount INT DEFAULT 0;
  DECLARE localResultRows INT DEFAULT 0;
  DECLARE localWorkId INT;
  DECLARE no_rows_indicator INT DEFAULT 0;
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char),', ',convert(IFNULL(inLowPosition,'NULL'),char),
  -- DEBUG_DIAGTABLE 	', ',convert(IFNULL(inHighPosition,'NULL'),char)));

  Create temporary table if not exists CachedWork (
	WorkId int auto_increment,
	PedPosId int, PedigreeSId varchar(32), PedTraitPosCM real,
	LC1MPId int, LC2MPId int, LC3MPId int,
	ServerId int, MarkerCount int,
	PRIMARY KEY (WorkId));

  set realLocusListType = inLocusListType;
  IF inLocusListType > 200 THEN
    set @MCMC_flag=1;
    set realLocusListType = (inLocusListType - 200) % 10;
    set totalSampleCount = (inLocusListType - 200) / 10; 
  ELSE
    set @MCMC_flag=0;
  END IF;
  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetWork: 1st try global @MCMC_flag is ', 
  -- DEBUG_DIAGTABLE   convert(IFNULL(@MCMC_flag,'NULL'),char),
  -- DEBUG_DIAGTABLE   ', local realLocusListType is ', convert(IFNULL(realLocusListType,'NULL'),char),
  -- DEBUG_DIAGTABLE   ', local totalSampleCount is ', convert(IFNULL(totalSampleCount,'NULL'),char)
  -- DEBUG_DIAGTABLE   ));

-- See if our cached work table has a row for us to return...

  SET no_rows_indicator = 0;
  Select WorkId, PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId
  into localWorkId, outPedPosId, outPedigreeSId, outPedTraitPosCM, outLC1MPId, outLC2MPId, outLC3MPId
  from CachedWork order by PedigreeSId limit 1;
  -- Crazy people have overloaded outLC2MPId to represent MarkerSetId. This is where bugs are born.
  IF realLocusListType = 1 THEN
    set @outMarkerSetId = outLC2MPId;
  END IF;
  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('Select from CachedWork got no_rows_indicator of ', convert(IFNULL(no_rows_indicator,'NULL'),char),
  -- DEBUG_DIAGTABLE   ', global @outMarkerSetId is ', convert(IFNULL(@outMarkerSetId,'NULL'),char)
  -- DEBUG_DIAGTABLE   ));

  IF no_rows_indicator THEN

    -- Someday we may want to verify that we didn't accidently change positions here and orphan a lot of work

    IF realLocusListType = 1 THEN
      call CacheMarkerWork(inServerId, inLowPosition, inHighPosition, localResultRows);
      set @outLC1MPId = -1;
    ELSE
      IF realLocusListType = 2 THEN
        call CacheTraitWork(inServerId, localResultRows);
      ELSE
        IF realLocusListType = 3 THEN
          call CacheCombinedWork(inServerId, inLowPosition, inHighPosition, localResultRows);
        ELSE
          IF realLocusListType = 100 THEN
            call Cache2ptInitialWork(inServerId, inLowPosition, inHighPosition, localResultRows);
          ELSE
            IF realLocusListType = 101 THEN
              call Cache2ptFinalWork(inServerId, inLowPosition, inHighPosition, localResultRows);
            -- DEBUG_DIAGTABLE ELSE
	    -- DEBUG_DIAGTABLE   Insert into Diag (Message) values (Concat('GetWork: unexpected inLocusListType of ',convert(IFNULL(inLocusListType,'NULL'),char)));
            END IF;
          END IF;
        END IF;
      END IF;
    END IF;
 
    -- Now try again...
    -- DEBUG_DIAGTABLE Insert into Diag (Message) values ('GetWork: 2nd attempt (after refreshing CachedWork');

    SET no_rows_indicator = 0;
    Select WorkId, PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId
    into localWorkId, outPedPosId, outPedigreeSId, outPedTraitPosCM, outLC1MPId, outLC2MPId, outLC3MPId
    from CachedWork order by PedigreeSId limit 1;
    IF realLocusListType = 1 THEN
      set @outMarkerSetId = outLC2MPId;
    END IF;
  END IF;
--  END IF;
  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetWork: 2nd try global @MCMC_flag is ', 
  -- DEBUG_DIAGTABLE   convert(IFNULL(@MCMC_flag,'NULL'),char),
  -- DEBUG_DIAGTABLE   ', local realLocusListType is ', convert(IFNULL(realLocusListType,'NULL'),char),
  -- DEBUG_DIAGTABLE   ', local totalSampleCount is ', convert(IFNULL(totalSampleCount,'NULL'),char)
  -- DEBUG_DIAGTABLE   ));

  IF no_rows_indicator THEN
    -- DEBUG_DIAGTABLE Insert into Diag (Message) values ('GetWork: no work found at all');
    Truncate table CachedWork;
  ELSE
    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetWork: returning ', convert(IFNULL(outPedPosId,'NULL'),char),', ', outPedigreeSId, ', ',
    -- DEBUG_DIAGTABLE convert(IFNULL(outPedTraitPosCM,'NULL'),char),', ', convert(IFNULL(outLC1MPId,'NULL'),char),', ', convert(IFNULL(outLC2MPId,'NULL'),char),', ', convert(IFNULL(outLC3MPId,'NULL'),char)));
    Delete from CachedWork where WorkId = localWorkId;

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

  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE localKeepAliveFlag INT;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('Cache2ptInitialWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char),', ',convert(IFNULL(inLowPosition,'NULL'),char),
  -- DEBUG_DIAGTABLE 	', ',convert(IFNULL(inHighPosition,'NULL'),char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  Cache2ptInitialWork: LOOP

    Select KeepAliveFlag into localKeepAliveFlag from Servers where ServerId = inServerId;
    IF localKeepAliveFlag = 0 THEN
      LEAVE Cache2ptInitialWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select fresh (no ServerId or StartTime), work which is generally to the right of the current marker
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount)
      select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId IS NULL
        limit 50 FOR UPDATE SKIP LOCKED; -- No ordering for best performance, SKIP LOCKED ignores already-locked rows
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values ('Cache2ptInitialWork: no work found at all');
      LEAVE Cache2ptInitialWork;
    ELSE
      -- Work we found has not been done at all, update the rows
      Select count(*) from CachedWork into outResultRows;
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('Cache2ptInitialWork: found ', convert(IFNULL(outResultRows,'NULL'),char), ' rows of undone work, marking'));
      Update Models a, CachedWork b set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId IS NULL;
      LEAVE Cache2ptInitialWork;
    END IF;

  END LOOP Cache2ptInitialWork;

  Commit;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('Cache2ptInitialWork: returning ', convert(IFNULL(outResultRows,'NULL'),char)));

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

  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE localKeepAliveFlag INT;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('Cache2ptFinalWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char),', ',convert(IFNULL(inLowPosition,'NULL'),char),
  -- DEBUG_DIAGTABLE   ', ',convert(IFNULL(inHighPosition,'NULL'),char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  Cache2ptFinalWork: LOOP

    Select KeepAliveFlag into localKeepAliveFlag from Servers where ServerId = inServerId;
    IF localKeepAliveFlag = 0 THEN
      LEAVE Cache2ptFinalWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select started but incomplete (has our ServerId, probably an entry in TP2MP) generally to the left of the current marker
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount)
      select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount
      from
	PedigreePositions PP, Models M, Servers S, TP2MP T
      where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.MarkerCount <= S.MarkerCount AND
	PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId = S.ServerId AND
	M.ModelId = T.ModelId AND
	M.EndTime IS NULL
        limit 50 FOR UPDATE SKIP LOCKED; -- No ordering for best performance, SKIP LOCKED ignores already-locked rows
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    LEAVE Cache2ptFinalWork;

  END LOOP Cache2ptFinalWork;

  Commit;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('Cache2ptFinalWork: returning ', convert(IFNULL(outResultRows,'NULL'),char)));

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

  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE localKeepAliveFlag INT;
  DECLARE localCandidatePedPosId INT;
  DECLARE localCandidateLimit INT;
  DECLARE localCandidateSMRT INT;
  DECLARE localCandidateFreeModels INT;
  DECLARE localExclusionCount INT;

  SET @dSString = "Not a SQL statement"; -- This cannot be a local DECLARE as prepared statement would be lost on exit, so I'm just flagging it.

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheCombinedWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char),', ',convert(IFNULL(inLowPosition,'NULL'),char),
  -- DEBUG_DIAGTABLE 	', ',convert(IFNULL(inHighPosition,'NULL'),char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;

  Truncate table CachedWork;

  Create temporary table if not exists ExcludedPedPosIds (PedPosId int);

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite, so be careful!
  CacheCombinedWork: LOOP

    Select KeepAliveFlag into localKeepAliveFlag from Servers where ServerId = inServerId;
    IF localKeepAliveFlag = 0 THEN
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values ('CacheCombinedWork: no longer kept alive');
      LEAVE CacheCombinedWork;
    END IF;

    SET outResultRows = 0;
    SET localCandidatePedPosId = NULL;
    SET localCandidateLimit = 51;

    -- The transaction can start AFTER the selection of a pedigree position to work on.

    -- Pick a random candidate PedigreePosition and limit Models to that PedPosId so we get a tolerable total runtime.
    -- Note that a server will have done all Trait and Marker Set models FIRST for any position, so we don't have to
    -- exclude them from our work selection.

    Select PP.PedPosId, PP.SingleModelRuntime, PP.FreeModels
      into localCandidatePedPosId, localCandidateSMRT, localCandidateFreeModels
   	from PedigreePositions PP, Servers S where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.PedPosId NOT IN (Select PedPosId from ExcludedPedPosIds) AND
	PP.FreeModels > 0 AND
	PP.MarkerCount = S.MarkerCount AND
	PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition
--        order by RAND() 
        limit 1; -- REMOVE for update

    IF localCandidatePedPosId IS NULL THEN
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values ('CacheCombinedWork: no work found');
      SET outResultRows = 0;
      LEAVE CacheCombinedWork;
    -- DEBUG_DIAGTABLE ELSE
    -- DEBUG_DIAGTABLE   Insert into Diag (Message) values (Concat('CacheCombinedWork: work with SMRT of ', convert(IFNULL(localCandidateSMRT,'NULL'),char), ' found for PedPosId ', convert(IFNULL(localCandidatePedPosId,'NULL'),char)));
    END IF;

    Start transaction;

    -- This is pointless for MC-MC since it'll do it 3000 times, so REALLY avoid the avg() call in Models
    IF @MCMC_flag = 1 THEN
      IF localCandidateSMRT IS NULL THEN
        Select avg(RunTimeCostSec) into localCandidateSMRT from Models where PedPosId=localCandidatePedPosId and RunTimeCostSec is not NULL;
      END IF;
      IF localCandidateSMRT IS NULL THEN
        SET localCandidateLimit = 5;
      ELSE
        IF localCandidateSMRT = 0 THEN
          SET localCandidateLimit = 51;
        ELSE
          SET localCandidateLimit = convert((60*60) / localCandidateSMRT, decimal(5,0));
          IF localCandidateLimit < 1 THEN
            SET localCandidateSMRT = 0;
            SET localCandidateLimit = 1;
          END IF;
          IF localCandidateLimit > 51 THEN
            SET localCandidateLimit = 51;
          END IF;
        END IF;
      END IF;
    ELSE
      SET localCandidateLimit = 51;
    END IF;
    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheCombinedWork: work for PedPosId ', convert(IFNULL(localCandidatePedPosId,'NULL'),char), ' with SMRT ', convert(IFNULL(localCandidateSMRT,'NULL'),char), ' limited to ', convert(IFNULL(localCandidateLimit,'NULL'),char)));
    Update Servers set CurrentPedPosId = localCandidatePedPosId, CurrentLimit = localCandidateLimit  where ServerId = inServerId;

    SET @dSString = Concat(
    	'Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount) ',
	'select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount ',
	'from PedigreePositions PP, Models M, Servers S ',
	'where S.ServerId = ', convert(IFNULL(inServerId,'NULL'),char),
	' AND PP.PedPosId = ', convert(IFNULL(localCandidatePedPosId,'NULL'),char),
	' AND M.PedPosId = ', convert(IFNULL(localCandidatePedPosId,'NULL'),char),
	' AND M.ServerId IS NULL limit ',
	convert(IFNULL(localCandidateLimit,'NULL'),char), ' FOR UPDATE SKIP LOCKED;');
    PREPARE dSHandle from @dSString;
    EXECUTE dSHandle;
    DEALLOCATE PREPARE dSHandle;

    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN
      Insert into ExcludedPedPosIds (PedPosId) value (localCandidatePedPosId);
      Select count(*) from ExcludedPedPosIds into localExclusionCount;
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheCombinedWork: no actual work (not ', convert(IFNULL(localCandidateFreeModels,'NULL'),char), ') found for ServerId ', convert(IFNULL(inServerId,'NULL'),char), ', PedPosId ', convert(IFNULL(localCandidatePedPosId,'NULL'),char), ', excluded ', convert(IFNULL(localExclusionCount,'NULL'),char), ', looping'));
      -- We don't leave, we try again with a different PedPosId (hopefully)
      Commit;
    ELSE
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheCombinedWork: found ', convert(IFNULL(outResultRows,'NULL'),char), ' rows of undone work, marking'));
      Update Models a, CachedWork b set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
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

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheCombinedWork: returning ', convert(IFNULL(outResultRows,'NULL'),char)));

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

  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE localKeepAliveFlag INT;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheMarkerWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char),', ',convert(IFNULL(inLowPosition,'NULL'),char),
  -- DEBUG_DIAGTABLE 	', ',convert(IFNULL(inHighPosition,'NULL'),char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  CacheMarkerWork: LOOP

    Select KeepAliveFlag into localKeepAliveFlag from Servers where ServerId = inServerId;
    IF localKeepAliveFlag = 0 THEN
      LEAVE CacheMarkerWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select first undone (no ServerId), and then if need be, less quality (lesser MarkerCount) work.
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount)
        select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, -1, MarkerSetId, -1, S.ServerId, S.MarkerCount
        from
	  PedigreePositions PP, MarkerSetLikelihood M, Servers S
        where
	  S.ServerId = inServerId AND
	  PP.StudyId = S.StudyId AND
	  PP.FreeModels > 0 AND
	  PP.MarkerCount <= S.MarkerCount AND
	  PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	  PP.ChromosomeNo = S.ChromosomeNo AND
	  PP.PedTraitPosCM > inLowPosition AND
	  PP.PedTraitPosCM <= inHighPosition AND
	  PP.PedPosId = M.PedPosId AND
	  M.ServerId IS NULL
	  limit 50 FOR UPDATE SKIP LOCKED; -- No ordering for best performance, SKIP LOCKED ignores already-locked rows
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN

      -- DEBUG_DIAGTABLE Insert into Diag (Message) values ('CacheMarkerWork: found no undone work, checking lesser quality work');

      SET outResultRows = 0;
        Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount)
          select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, -1, MarkerSetId, -1, M1.ServerId, S.MarkerCount
          from
	    PedigreePositions PP, MarkerSetLikelihood M1, Servers S
          where
	    S.ServerId = inServerId AND
	    PP.StudyId = S.StudyId AND
	    PP.PedPosId = M1.PedPosId AND
	    PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	    PP.ChromosomeNo = S.ChromosomeNo AND
	    PP.PedTraitPosCM > inLowPosition AND
	    PP.PedTraitPosCM <= inHighPosition AND
	    PP.MarkerCount <= S.MarkerCount AND
	    M1.MarkerCount < S.MarkerCount AND
	    NOT EXISTS (
		Select * from MarkerSetLikelihood M2 where 
		M1.PedPosId = M2.PedPosId AND
		M2.MarkerCount = S.MarkerCount
		)
	    limit 50 FOR UPDATE SKIP LOCKED; -- No ordering for best performance, SKIP LOCKED ignores already-locked rows
      Select count(*) from CachedWork into outResultRows;

      IF outResultRows = 0 THEN
	-- DEBUG_DIAGTABLE Insert into Diag (Message) values ('CacheMarkerWork: no work found at all');
        LEAVE CacheMarkerWork;
      ELSE
	-- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheMarkerWork: found ', convert(IFNULL(outResultRows,'NULL'),char), ' rows of lesser-quality, prepping by inserting new Models row'));
        Insert into MarkerSetLikelihood (PedPosId, ServerId, MarkerCount, StartTime)
	 select PedPosId, ServerId, MarkerCount, CURRENT_TIMESTAMP from CachedWork;
	select LAST_INSERT_ID() into @outMarkerSetId;
        LEAVE CacheMarkerWork;
      END IF;

    ELSE
      -- Work we found has not been done at all, update the rows
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheMarkerWork: found ', convert(IFNULL(outResultRows,'NULL'),char), ' rows of undone work, marking'));
      Select count(*) from CachedWork into outResultRows;
      Update MarkerSetLikelihood a, CachedWork b  set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
       	a.PedPosId = b.PedPosId AND
	b.LC1MPId = -1 AND
	b.LC3MPId = -1 AND
	a.ServerId IS NULL;
      LEAVE CacheMarkerWork;
    END IF;

  END LOOP CacheMarkerWork;

  Commit;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheMarkerWork: returning ', convert(IFNULL(outResultRows,'NULL'),char)));

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

  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE localKeepAliveFlag INT;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheTraitWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char)));

  Update Servers set LastHeartbeat = CURRENT_TIMESTAMP where ServerId = inServerId;
  Commit;

  Truncate table CachedWork;

  -- This is a dummy loop that allows us to exit whenever we want. Its only a dummy as long as all paths do exit,
  -- otherwise, it's infinite!
  CacheTraitWork: LOOP

    Select KeepAliveFlag into localKeepAliveFlag from Servers where ServerId = inServerId;
    IF localKeepAliveFlag = 0 THEN
      LEAVE CacheTraitWork;
    END IF;
    Commit;

    Start transaction;

    SET outResultRows = 0;
    -- Select first undone (no ServerId), and then if need be, less quality (lesser MarkerCount) work.
    Insert into CachedWork (PedPosId, PedigreeSId, PedTraitPosCM, LC1MPId, LC2MPId, LC3MPId, ServerId, MarkerCount)
      select PP.PedPosId, PP.PedigreeSId, PP.PedTraitPosCM, M.LC1MPId, M.LC2MPId, M.LC3MPId, S.ServerId, S.MarkerCount
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND
	PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.RefTraitPosCM = -9999.99 AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId IS NULL
        limit 50 FOR UPDATE SKIP LOCKED; -- No ordering for best performance, SKIP LOCKED ignores already-locked rows
--      order by ModelId limit 10; -- Test version that keeps the models in-order
--      order by RAND() limit 10; -- The RAND() bit keeps identical servers from fighting too much.
    Select count(*) from CachedWork into outResultRows;

    IF outResultRows = 0 THEN
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values ('CacheTraitWork: found no undone work (none at all for traits!)');
      LEAVE CacheTraitWork;
    ELSE
      -- Work we found has not been done at all, update the rows
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheTraitWork: found ', convert(IFNULL(outResultRows,'NULL'),char), ' rows of undone work, marking'));
      Select count(*) from CachedWork into outResultRows;
      Update Models a, CachedWork b set a.ServerId = b.ServerId, a.MarkerCount = b.MarkerCount, a.StartTime = CURRENT_TIMESTAMP where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId IS NULL;
      LEAVE CacheTraitWork;
    END IF;

  END LOOP CacheTraitWork;

  Commit;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CacheTraitWork: returning ', convert(IFNULL(outResultRows,'NULL'),char)));

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

  DECLARE markerWork INT DEFAULT 0;
  DECLARE totalWork INT DEFAULT 0;
  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE no_rows_indicator INT DEFAULT 0;
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CountWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char),', ',convert(IFNULL(inLowPosition,'NULL'),char),
  -- DEBUG_DIAGTABLE 	', ',convert(IFNULL(inHighPosition,'NULL'),char)));

  -- make sure we have trait position in the given range first before the expensive join
  SET no_rows_indicator = 0;
  select count(*) into outWorkCount from PedigreePositions PP, Servers S
  where 
	S.ServerId = inServerId AND
        PP.StudyId = S.StudyId AND
	PP.MarkerCount <= S.MarkerCount AND
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition;

  IF outWorkCount != 0 THEN
  -- Select first undone (no ServerId), and then if need be, less quality (lesser MarkerCount) work.

  SET no_rows_indicator = 0;
  SET outWorkCount = 0;
  Select count(*) into outWorkCount
  from
      PedigreePositions PP, Models M, Servers S
  where
	S.ServerId = inServerId AND
        PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId IS NULL;
    Select count(*) into markerWork
    from
        PedigreePositions PP, MarkerSetLikelihood M, Servers S
    where 
	S.ServerId = inServerId AND
        PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M.PedPosId AND
	M.ServerId IS NULL;
  SET outWorkCount = outWorkCount + markerWork;
  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CountWork: found ', convert(IFNULL(outWorkCount,'NULL'),char),' undone work'));
  IF outWorkCount= 0 THEN
    Select count(*) into outWorkCount
    from
        PedigreePositions PP, Models M1, Servers S
    where
        S.ServerId = inServerId AND
        PP.StudyId = S.StudyId AND
	PP.MarkerCount <= S.MarkerCount AND
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
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
    Select count(*) into markerWork
    from	
        PedigreePositions PP, MarkerSetLikelihood M1, Servers S
    where 
	S.ServerId = inServerId AND
        PP.StudyId = S.StudyId AND
	PP.FreeModels > 0 AND
	PP.MarkerCount <= S.MarkerCount AND
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
	PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedTraitPosCM > inLowPosition AND
	PP.PedTraitPosCM <= inHighPosition AND
	PP.PedPosId = M1.PedPosId AND
	M1.MarkerCount < S.MarkerCount AND
	NOT EXISTS (
		Select * from Models M2 where 
		M1.PedPosId = M2.PedPosId AND
		M2.MarkerCount = S.MarkerCount
		);

    SET outWorkCount = outWorkCount + markerWork;
    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CountWork: found ', convert(IFNULL(outWorkCount,'NULL'),char),' lower-quality  work'));
  END IF;
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

  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE inSampleId int;
  DECLARE sampleCount int;
  DECLARE localOldSMRT INT;
  DECLARE localNewSMRT INT;
  DECLARE localOldLimit INT;
  DECLARE localOutModelId INT;
  DECLARE localOutWeightedLRComponent DOUBLE;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('PutWork: called w/ ', convert(IFNULL(inServerId,'NULL'),char), ', ', convert(IFNULL(inPedPosId,'NULL'),char),', ',
  -- DEBUG_DIAGTABLE 	convert(IFNULL(inLC1MPId,'NULL'),char),', ', convert(IFNULL(inLC2MPId,'NULL'),char),', ', convert(IFNULL(inLC3MPId,'NULL'),char), ', ',
  -- DEBUG_DIAGTABLE 	convert(IFNULL(inMarkerCount,'NULL'),char),', ', convert(IFNULL(inLikelihood,'NULL'),char),', ', convert(IFNULL(inRuntimeCostSec,'NULL'),char)));

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('PutWork: @outMarkerSetId ', convert(IFNULL(@outMarkerSetId,'NULL'), char)));

  Set localNewSMRT = 0;

  -- This doesn't have to be in the transaction
  Select ModelId INTO localOutModelId from Models where
	PedPosId = inPedPosId AND
	LC1MPId = inLC1MPId AND
	LC2MPId = inLC2MPId AND
	LC3MPId = inLC3MPId AND
	ServerId = inServerId;

  start transaction; 

  -- NO LONGER lock the PedigreePosition row first to avoid deadlock with CacheCombinedWork.
  -- Both routines have the same locking order this way, with PedigreePositions table first.
  -- Updating likelihood should HAVE beEN quick. BUT MCMC DOES A GAZILLION UPDATES.

  -- Select PedPosId into inPedPosId from PedigreePositions where PedPosId = inPedPosId for update;

  IF inMarkerCount = 100 THEN
    -- 2pt, only have the initial first half, hold onto it.
    Insert into TP2MP (ModelId, WeightedLRComponent, RuntimeCostSec) values (localOutModelId, inLikelihood, inRuntimeCostSec);
  ELSE
    IF inMarkerCount = 101 THEN
      -- 2pt, got the final half, combine them and truely finish the job.
      Select WeightedLRComponent, RuntimeCostSec into localOutWeightedLRComponent, localNewSMRT from TP2MP where ModelId = localOutModelId;
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('PutWork: Adding inLikelihood ', convert(IFNULL(inLikelihood,'NULL'),char), ' to existing ', convert(IFNULL(localOutWeightedLRComponent,'NULL'),char)));
      Update Models set Likelihood = inLikelihood+localOutWeightedLRComponent, RuntimeCostSec = inRuntimeCostSec+localNewSMRT,
	MarkerCount = inMarkerCount, EndTime = CURRENT_TIMESTAMP where ModelId = localOutModelId;
    ELSE
      -- Classic multipoint
      IF @outLC1MPId >= 0 THEN
        -- trait or combined likelihood
        Update Models set Likelihood = inLikelihood, RuntimeCostSec = inRuntimeCostSec, MarkerCount = inMarkerCount, EndTime = CURRENT_TIMESTAMP where
	  ModelId = localOutModelId;
      ELSE
        -- Markerset likelihood
	IF @MCMC_flag = 0 THEN
	  Update MarkerSetLikelihood set Likelihood = inLikelihood, RuntimeCostSec = inRuntimeCostSec, MarkerCount = inMarkerCount, EndTime = CURRENT_TIMESTAMP where
	  PedPosId = inPedPosId AND ServerId = inServerId;
	ELSE
	  -- MCMC weirdness
	  set inSampleId = (inMarkerCount - 200) / 10; 
	  IF inSampleId > 0 THEN
            select count(*) into sampleCount from MarkerSetLikelihood_MCMC where MarkerSetId=@outMarkerSetId and SampleId = inSampleId;
	    IF sampleCount = 0 THEN
              -- insert
              insert into MarkerSetLikelihood_MCMC (MarkerSetId, SampleId, Likelihood) values (@outMarkerSetId, inSampleId, inLikelihood);
	    ELSE
	      -- update
	      Update MarkerSetLikelihood_MCMC set Likelihood = inLikelihood where MarkerSetId=@outMarkerSetId and SampleId = inSampleId;
            END IF;
	  ELSE
            -- set the marker likelihood to 1
	    Update MarkerSetLikelihood set Likelihood = 1, RunTimeCostSec = inRuntimeCostSec, MarkerCount = inMarkerCount, EndTime = CURRENT_TIMESTAMP where
	    PedPosId = inPedPosId AND ServerId = inServerId;
          END IF;
	END IF;
      END IF;
    END IF;
  -- Clear-out any RegionModels so we know they're not pending
  -- 
  -- Note that we're not actually regularly using RegionModels nowadays, but it
  -- can get used now and again in debug situations, and, well, better safe
  -- than sorry.
  Delete from RegionModels where ModelId = localOutModelId;
  END IF;
--  END IF;

  -- Now the hard part. If this model took an unexpectedly long time, we need to update the PedPosId with that fact, and
  -- we want to jettison the remaining models because  we don't want to be holding lots of work.
  -- If you wonder why we have to do this, look at the variation on SingleModelRuntimes using a query like this:
  -- Select PedPosId, count(*), max(RuntimeCostSec), min(RuntimeCostSec), avg(RuntimeCostSec), stddev_pop(RuntimeCostSec) from Models where PedPosId in (Select PedPosId from PedigreePositions where StudyId = 8) group by PedPosId order by max(RuntimeCostSec) desc limit 20;
  IF inRuntimeCostSec > 60 THEN
    Select SingleModelRuntime into localOldSMRT from PedigreePositions where PedPosId = inPedPosId;
    IF localOldSMRT IS NULL THEN
      SET localOldSMRT = 0;
    END IF;
    IF (inRuntimeCostSec * 0.70) > localOldSMRT THEN
      -- Set the SingleModelRuntime to the new value
      Update PedigreePositions set SingleModelRuntime = inRuntimeCostSec where PedPosId = inPedPosId;
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('PutWork: Reset SMRT due to ', convert(IFNULL(inRuntimeCostSec,'NULL'),char), 's model for PedPosId ', convert(IFNULL(inPedPosId,'NULL'),char), ', which had ', convert(IFNULL(localOldSMRT,'NULL'),char), 's SMRT'));
    END IF;
    -- We can't just look at the oldSMRT because we could be holding work for which some other server just updated the cost.
    Select inRuntimeCostSec * count(*) into localOldLimit from CachedWork;
    IF localOldLimit > (60*60*0.70) THEN
      -- Release the models we haven't processed
      Update Models a, CachedWork b  set a.ServerId = NULL, a.MarkerCount = b.MarkerCount, a.StartTime = NULL where
	a.PedPosId = b.PedPosId AND
	a.LC1MPId = b.LC1MPId AND
	a.LC2MPId = b.LC2MPId AND
	a.LC3MPId = b.LC3MPId AND
	a.ServerId = inServerId;
      -- Clear-out everything so we get a new, fresh bunch of PROPERLY SIZED work
      Delete from CachedWork;
      -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('PutWork: flushed remaining work due to ', convert(IFNULL(inRuntimeCostSec,'NULL'),char), 's model for PedPosId ', convert(IFNULL(inPedPosId,'NULL'),char), ', which had ', convert(IFNULL(localOldSMRT,'NULL'),char), 's SMRT'));
    END IF;
  END IF;

  commit;
END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS SetDummyNullLikelihood;
DELIMITER //

CREATE PROCEDURE SetDummyNullLikelihood (IN inServerId int)
BEGIN

  DECLARE version char(96) DEFAULT '$Id$';

-- Set all trait and marker likelhood Models for pedigrees handled by this ServerId to be finished and have an Likelihood of 1.0

  DECLARE inModelId int;
  DECLARE inPedPosId int;
  DECLARE EmptyCursor BOOLEAN DEFAULT 0;
  DECLARE DummyMarkerModels CURSOR FOR
      Select M.PedPosId
      from
	PedigreePositions PP, MarkerSetLikelihood M, Servers S
      where
	S.ServerId = inServerId AND PP.StudyId = S.StudyId AND PP.FreeModels > 0 AND PP.MarkerCount = 1 AND -- WANT TO NOT CONSIDER FREE MODELS ANYMORE...
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
        PP.ChromosomeNo = S.ChromosomeNo AND
	PP.PedPosId = M.PedPosId AND M.ServerId IS NULL;
  DECLARE DummyTraitModels CURSOR FOR
      Select M.ModelId
      from
	PedigreePositions PP, Models M, Servers S
      where
	S.ServerId = inServerId AND PP.StudyId = S.StudyId AND PP.FreeModels > 0 AND PP.MarkerCount = 1 AND -- WANT TO NOT CONSIDER FREE MODELS ANYMORE...
        PP.PedigreeSId in (select PedigreeSId from ServerPedigrees where ServerId=inServerId) and
        PP.ChromosomeNo = S.ChromosomeNo AND
	PP.RefTraitPosCM = -9999.99 AND PP.PedPosId = M.PedPosId AND M.ServerId IS NULL;
  DECLARE CONTINUE HANDLER FOR SQLSTATE '02000' SET EmptyCursor = 1;

  -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('SetDummyNullLikelihood: called w/ ', convert(IFNULL(inServerId,'NULL'),char)));

  OPEN DummyMarkerModels;
  REPEAT
    FETCH DummyMarkerModels into inPedPosId;
    Update MarkerSetLikelihood set ServerId = inServerId, StartTime = CURRENT_TIMESTAMP, Likelihood = 1.0, EndTime = CURRENT_TIMESTAMP where PedPosId=inPedPosId;
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

    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('GetDLikelihood: returning.'));

END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS CleanOrphans;
DELIMITER //

CREATE PROCEDURE CleanOrphans (IN inStudyLabel varchar(64))
BEGIN

  DECLARE inStudyId int;
  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE no_rows_indicator INT DEFAULT 0;

  -- Put orphaned Models back up for adoption.

  DECLARE inModelId int;
  DECLARE inPedPosId int;
  DECLARE inMarkerCount int;
  DECLARE orphanCount int DEFAULT 0;

  DECLARE EmptyCursor BOOLEAN DEFAULT 0;
  DECLARE OrphanModels CURSOR FOR
    Select ModelId, PedPosId from Models where
	StartTime is NOT NULL AND Likelihood is NULL AND
	PedPosId in (Select PedPosId from PedigreePositions where StudyId = inStudyId) AND
	ServerId NOT IN (Select ServerId from Servers where ExitStatus IS NULL);
  DECLARE OrphanMarkers CURSOR FOR
    Select PedPosId, MarkerCount from MarkerSetLikelihood where
	StartTime is NOT NULL AND Likelihood is NULL AND
	PedPosId in (Select PedPosId from PedigreePositions where StudyId = inStudyId) AND
	ServerId NOT IN (Select ServerId from Servers where ExitStatus IS NULL);
  DECLARE CONTINUE HANDLER FOR SQLSTATE '02000' SET EmptyCursor = 1;

-- Get the actual StudyId
  WholeThing: LOOP
  SET no_rows_indicator = 0;
  Select StudyId into inStudyId from Studies where StudyLabel = inStudyLabel;
  IF no_rows_indicator THEN
-- DEBUG_DIAGTABLE     Insert into Diag (Message) values (Concat('CleanOrphans(', inStudyLabel, '): cannot find the study!'));
      Leave WholeThing;      
  END IF;


  OPEN OrphanModels;
  FETCH OrphanModels into inModelId, inPedPosId;
  WHILE EmptyCursor <> 1 DO

    UPDATE Models set ServerId = NULL, StartTime = NULL where ModelId = inModelId;
    DELETE from TP2MP where ModelId = inModelId;
    SET orphanCount = orphanCount+1;
    FETCH OrphanModels into inModelId, inPedPosId;

  END WHILE;
  CLOSE OrphanModels;

  set EmptyCursor=0;
  OPEN OrphanMarkers;
  FETCH OrphanMarkers into inPedPosId, inMarkerCount;
  WHILE EmptyCursor <> 1 DO

    UPDATE MarkerSetLikelihood set ServerId = NULL, StartTime = NULL where PedPosId=inPedPosId and MarkerCount=inMarkerCount;
    SET orphanCount = orphanCount+1;
    FETCH OrphanMarkers into inPedPosId, inMarkerCount;

  END WHILE;
  CLOSE OrphanMarkers;

  Select orphanCount 'Orphans Freed';
  Leave WholeThing;
  END LOOP WholeThing; 

END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS CleanStudy;
DELIMITER //

CREATE PROCEDURE CleanStudy (IN inStudyLabel varchar(64))
BEGIN

  DECLARE inStudyId int;
  DECLARE version char(96) DEFAULT '$Id$';
  DECLARE no_rows_indicator INT DEFAULT 0;
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

  WholeThing: LOOP

-- Get the actual StudyId
  SET no_rows_indicator = 0;
  Select StudyId into inStudyId from Studies where StudyLabel = inStudyLabel;
  IF no_rows_indicator THEN
    -- DEBUG_DIAGTABLE Insert into Diag (Message) values (Concat('CleanStudy(', inStudyLabel, '): cannot find the study!'));
    LEAVE WholeThing;
  END IF;

-- Remove all rows associated with a study

  Delete from ServerPedigrees where ServerId in (Select ServerId from Servers where StudyId = inStudyId);
  Delete from MarkerSetLikelihood_MCMC where MarkerSetId in (Select MarkerSetId from MarkerSetLikelihood where PedPosId in (Select PedPosId from PedigreePositions where StudyId = inStudyId));
  Delete from MarkerSetLikelihood where PedPosId in (Select PedPosId from PedigreePositions where StudyId = inStudyId);
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
  Delete from Analyses where StudyId = inStudyId;
  Delete from Studies where StudyId = inStudyId;
  -- These next three statements are for tables we don't currently use, but
  -- better safe than sorry
  Delete from RegionModels where RegionId in (Select distinct a.RegionId from Regions a, Analyses b where a.AnalysisId = b.AnalysisId AND b.StudyId = inStudyId);
  Delete from SingleSizingRuns where StudyId = inStudyId;
  Delete from SingleModelRuntimes where StudyId = inStudyId;

  Leave WholeThing;
	
END LOOP WholeThing;

END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS MoveStudy;
DELIMITER //

CREATE PROCEDURE MoveStudy (IN inFromStudyId int, IN inToStudyId int)
BEGIN

  DECLARE version char(96) DEFAULT '$Id$';

-- Move all rows associated with a StudyId to a new StudyId without violating integrity constraints
-- NEW STUDYID MUST NOT ALREADY EXIST

  Update Studies set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update Analyses set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update Maps set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update Markers set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update MapMarkers set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update Pedigrees set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update Positions set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update PedigreePositions set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update Regions set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update Servers set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update LGModels set StudyId = inToStudyId where StudyId = inFromStudyId;
  -- These next two statements are for tables we don't currently use, but
  -- better safe than sorry
  Update SingleSizingRuns set StudyId = inToStudyId where StudyId = inFromStudyId;
  Update SingleModelRuntimes set StudyId = inToStudyId where StudyId = inFromStudyId;

END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS UnloadModelParts;
DELIMITER //

CREATE PROCEDURE UnloadModelParts ()
BEGIN

  DECLARE version char(96) DEFAULT '$Id$';

  Select Concat(
    'Insert into DModelParts (MPId, DGF, BigPen, BigLittlePen, LittleBigPen, LittlePen) values (',
	convert(IFNULL(MPId,'NULL'),char),',',
	format(DGF,32),',',
	format(BigPen,32),',',
	format(BigLittlePen,32),',',
	format(LittleBigPen,32),',',
	format(LittlePen,32),');')
    from DModelParts;

  Select Concat(
    'Insert into QModelParts (MPId, DGF, BigMean, BigLittleMean, LittleBigMean, LittleMean, BigSD, BigLittleSD, LittleBigSD, LittleSD, Threshold) values (',
	convert(IFNULL(MPId,'NULL'),char),',',
	format(DGF,32),',',
	format(BigMean,32),',',
	format(BigLittleMean,32),',',
	format(LittleBigMean,32),',',
	format(LittleMean,32),',',
	format(BigSD,32),',',
	format(BigLittleSD,32),',',
	format(LittleBigSD,32),',',
	format(LittleSD,32),',',
	format(Threshold,32),');')
    from QModelParts;

END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS Reconcile;
DELIMITER //

CREATE PROCEDURE Reconcile(IN inAnalysisId int)
-- Reconcile: Mark any servers not really in processlist with ExitStatus 42 and
-- free up any reserved models for servers that failed for any reason
  -- This was pulled out of 'Q' so we could parameterize with the analysis ID,
  -- which is used to constrain the list of servers (otherwise, we get a lot of
  -- lockups, particularly when doing byped analysis). It's used by cycle
  -- scripts.
BEGIN
  -- Mark any servers that never got the chance to report doom
  -- UPDATE Servers a, Analyses b SET a.ExitStatus = 42 WHERE a.ConnectionId NOT IN (SELECT ID FROM INFORMATION_SCHEMA.PROCESSLIST) AND a.ExitStatus IS NULL AND a.PedigreeRegEx = b.PedigreeRegEx AND a.PedigreeNotRegEx = b.PedigreeNotRegEx AND b.AnalysisId = inAnalysisId;
  -- This use of a temporary table is a tad silly. We *would* just go with
  -- what's above rather than using a temporary table as a go-between, but
  -- apparently when you have a sufficiently large number of connections
  -- getting started at a given time MySQL goes all "omg wtf nooooooo" and
  -- gives us the singularly unhelpful error message:
  -- Data too long for column 'USER' at row 1
  -- A lot of trial and error revealed that putting in a temporary table
  -- for the subquery "fixes" it. Thanks, MySQL!
  DROP TABLE IF EXISTS silly_mysql_bug_workaround_for_reconcile;
  CREATE TEMPORARY TABLE silly_mysql_bug_workaround_for_reconcile AS SELECT ID FROM INFORMATION_SCHEMA.PROCESSLIST;
  UPDATE Servers a, Analyses b SET a.ExitStatus = 42 WHERE a.ConnectionId NOT IN (SELECT ID FROM silly_mysql_bug_workaround_for_reconcile) AND a.ExitStatus IS NULL AND a.PedigreeRegEx = b.PedigreeRegEx AND a.PedigreeNotRegEx = b.PedigreeNotRegEx AND b.AnalysisId = inAnalysisId;
  -- Free all still reserved models for which the server did not exit successfully.
  UPDATE Models a, Servers b, Analyses c SET a.ServerId = NULL, a.StartTime = NULL WHERE a.ServerId = b.ServerId AND b.ExitStatus IS NOT NULL AND a.Likelihood IS NULL AND b.PedigreeRegEx = c.PedigreeRegEx AND b.PedigreeNotRegEx = c.PedigreeNotRegEx AND c.AnalysisId = inAnalysisId;
END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS FixFree;
DELIMITER //

CREATE PROCEDURE FixFree(IN inAnalysisId int)
-- FixFree: Turn on or off FreeModels flags for appropriate PedigreePositions
  -- Pulled out of 'Q' so it can be parameterized with the analysis ID - needed
  -- to keep this from causing a lot of lockup issues when it does updates.
  -- It's used by cycle scripts.
BEGIN
  DROP TABLE IF EXISTS FreeModelFlags;
  CREATE TEMPORARY TABLE FreeModelFlags AS SELECT a.PedPosId, 0 StillFreeModels FROM
  PedigreePositions a, Analyses b WHERE a.PedigreeSId RLIKE b.PedigreeRegEx AND a.PedigreeSId NOT RLIKE b.PedigreeNotRegEx AND b.AnalysisId = inAnalysisId;
  UPDATE FreeModelFlags a, Models b SET a.StillFreeModels = 1 WHERE a.PedPosId = b.PedPosId AND b.Likelihood IS NULL;
  UPDATE PedigreePositions a, FreeModelFlags b SET a.FreeModels = b.StillFreeModels WHERE a.PedPosId = b.PedPosId;
END;
//
DELIMITER ;

DROP PROCEDURE IF EXISTS Q;
DELIMITER //

CREATE PROCEDURE Q (IN inWhich varchar(16))
-- A series of predefined helpful queries.
BEGIN

  DECLARE localTotal int;
  DECLARE localIgnore int;
  DECLARE localStartTime TIME;

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
      Select SubTime(Now(), '01:00') into localStartTime from dual;
      Select count(*) into localTotal from Servers s, Models m where m.ServerId = s.ServerId AND m.Likelihood is NOT NULL AND s.LastHeartbeat > localStartTime;
      Select Sleep (60) into localIgnore;
      Select Now(), (count(*) - localTotal)/60.0 'MPS' from Servers s, Models m where m.ServerId = s.ServerId AND m.Likelihood is NOT NULL AND s.LastHeartbeat > localStartTime;
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

  -- TotalFree: Show total unallocated work by StudyLabel
  IF inWhich = 'TotalFree' THEN
    Select S.StudyId, S.StudyLabel, P.PedigreeSId, count(*) from
      Studies S, Pedigrees P, PedigreePositions PP, Models M where
      S.StudyId = P.StudyId AND P.PedigreeSId = PP.PedigreeSId AND S.StudyId = PP.StudyId AND
      PP.PedPosId = M.PedPosId AND M.ServerId is NULL group by S.StudyLabel, P.PedigreeSId;
    Leave WholeThing;
  END IF;

  -- Reconcile: Mark any servers not really in processlist with ExitStatus 42
  -- and free up any reserved models for servers that failed for any reason
  -- DEPRECATED - this should ONLY be used if nothing else is interacting with
  -- the server, as it can otherwise cause serious contention issues because of
  -- the Update statements (see Reconcile(analysisid) above, instead)
  IF inWhich = 'Reconcile' THEN
    -- Mark any servers that never got the chance to report doom
    -- Update Servers set ExitStatus = 42 where ConnectionId NOT IN (Select ID from INFORMATION_SCHEMA.PROCESSLIST) AND ExitStatus IS NULL;
    -- This use of a temporary table is a tad silly. We *would* just go with
    -- what's above rather than using a temporary table as a go-between, but
    -- apparently when you have a sufficiently large number of connections
    -- getting started at a given time MySQL goes all "omg wtf nooooooo" and
    -- gives us the singularly unhelpful error message:
    -- Data too long for column 'USER' at row 1
    -- A lot of trial and error revealed that putting in a temporary table
    -- for the subquery "fixes" it. Thanks, MySQL!
    Create temporary table silly_mysql_bug_workaround_for_reconcile as Select ID from INFORMATION_SCHEMA.PROCESSLIST;
    Update Servers set ExitStatus = 42 where ConnectionId NOT IN (Select ID from silly_mysql_bug_workaround_for_reconcile) AND ExitStatus IS NULL;
    DROP TABLE IF EXISTS silly_mysql_bug_workaround_for_reconcile;
    -- Free all still reserved models for which the server did not exit
    -- successfully.
    Update Models a, Servers b set a.ServerId = NULL, a.StartTime = NULL where a.ServerId = b.ServerId AND b.ExitStatus IS NOT NULL AND a.Likelihood IS NULL;
 
    Leave WholeThing;
  END IF;
  -- FixFree: Turn on or off FreeModels flags for appropriate PedigreePositions
  -- DEPRECATED - this should ONLY be used if nothing else is interacting with
  -- the server, as it can otherwise cause serious contention issues because of
  -- the Update statements (see FixFree(analysisid) above, instead)
  IF inWhich = 'FixFree' THEN
    Update PedigreePositions set FreeModels = 0 where FreeModels <> 0;
    Create temporary table FreeModelCounts as Select PedPosId, count(*) 'FreeModels' from Models where Likelihood IS NULL group by PedPosId;
    Update PedigreePositions a, FreeModelCounts b set a.FreeModels = b.FreeModels where a.PedPosId = b.PedPosId;
    Drop table FreeModelCounts;
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
	('TotalFree', 'Show total unallocated work by StudyLabel'),
	('Reconcile', 'Mark any servers not really in processlist with ExitStatus 42 and release any incomplete work (DEPRECATED, use Reconcile(analysisid) instead)'),
	('FixFree', 'Reset FreeModels flags (DEPRECATED, use FixFree(analysisid) instead)');
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

DROP PROCEDURE IF EXISTS ServerSignOn;
DELIMITER //

CREATE PROCEDURE ServerSignOn (IN inHostName varchar(32), IN inProcessId int, IN inKeepAliveFlag int,
	IN inStudyId int, IN inPedigreeRegEx varchar(1025), IN inPedigreeNotRegEx varchar(1025), 
	IN inChromosomeNo int, IN inAlgorithm varchar(2), IN inMarkerCount int, 
	IN inProgramVersion varchar(32), IN inSampleIdStart int, IN inSampleIdEnd int, 
	OUT outServerId int)
BEGIN
  start transaction;
    Insert into Servers ( ConnectionId, HostName, ProcessId, KeepAliveFlag, StudyId, PedigreeRegEx, PedigreeNotRegEx, ChromosomeNo, Algorithm, MarkerCount, ProgramVersion, SampleIdStart, SampleIdEnd) values (connection_id(),inHostName, inProcessId, inKeepAliveFlag, inStudyId, inPedigreeRegEx, inPedigreeNotRegEx, inChromosomeNo, inAlgorithm, inMarkerCount, inProgramVersion, inSampleIdStart, inSampleIdEnd);
    select LAST_INSERT_ID() into outServerId;
    -- Figure out what pedigrees satisfy the regular expressions 
    -- using ServerPedigrees improves the query performance for GetWork as
    -- the locking will be on specific sets of pedigrees, instead of all the pedigrees
    -- if using RLIKE in those select for updates
    Insert into ServerPedigrees (ServerId, PedigreeSId) 
	select outServerId, PedigreeSId from Pedigrees where StudyId=inStudyId and 
          PedigreeSId RLIKE inPedigreeRegEx and 
          PedigreeSId NOT RLIKE inPedigreeNotRegEx;

  commit;
END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS ServerSignOff;
DELIMITER //

CREATE PROCEDURE ServerSignOff (IN inServerId int, IN inExitStatus int)
BEGIN
  start transaction;
    delete from ServerPedigrees where ServerId=inServerId;
    update Servers set ExitStatus=inExitStatus, StopTime=NOW() where ServerId=inServerId;
  commit;
END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetStudyId;
DELIMITER //

-- StudyLabel, LiabilityClassCnt and ImprintingFlag combination should be unique
CREATE PROCEDURE GetStudyId (IN inStudyLabel varchar(64), IN inLiabilityClassCnt int, IN inImprintingFlag char(1), OUT outStudyId int)
BEGIN
  DECLARE version char(96) DEFAULT '$Id$ in personal kelvin/trunk/LKS';
  DECLARE no_rows_indicator INT DEFAULT 0; 
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

  Select StudyId into outStudyId from Studies where StudyLabel=inStudyLabel and LiabilityClassCnt=inLiabilityClassCnt and ImprintingFlag=inImprintingFlag;
  IF no_rows_indicator THEN
    start transaction;
      insert into Studies (StudyLabel, LiabilityClassCnt, ImprintingFlag) values (inStudyLabel, inLiabilityClassCnt, inImprintingFlag);
      select LAST_INSERT_ID() into outStudyId;
    commit;
  END IF;
END
//
DELIMITER ;

DROP PROCEDURE IF EXISTS GetMapId;
DELIMITER //

-- inReferenceFlag 1-Reference/client map 0-StudySpecific/Server map
CREATE PROCEDURE GetMapId (IN inStudyId int, IN inReferenceFlag int, IN inMapScale char(1), IN inDescription varchar(128), OUT outMapId int, OUT outMapScale char(1), OUT outDescription varchar(128))
BEGIN
  DECLARE no_rows_indicator INT DEFAULT 0; 
  DECLARE no_rows CONDITION FOR 1329;
  DECLARE CONTINUE HANDLER FOR no_rows SET no_rows_indicator = 1;

  set outMapScale=" ";
  set outDescription="";
  -- dealing with a client/reference map
  -- see whether we already have a reference map set for this study
  IF inReferenceFlag=1 THEN
    -- client/ref map
    select m.MapId, m.MapScale, m.Description into outMapId, outMapScale, outDescription from Maps m, Studies s  where m.ReferenceFlag=1 and m.StudyId=inStudyId and m.StudyId=s.StudyId and s.ReferenceMapId=m.MapId;
    IF no_rows_indicator THEN
      -- this is the first time we add in a reference map for this study
      start transaction;
      insert into Maps (StudyId, MapScale, Description, ReferenceFlag) values (inStudyId, inMapScale, inDescription, inReferenceFlag);
      select LAST_INSERT_ID() into outMapId;
      -- update the study table as well
      update Studies set ReferenceMapId=outMapId where StudyId=inStudyId;
      commit;
    ELSE
      -- the reference map is already set, check whether the file match
      IF outDescription<>inDescription THEN
        -- conflict of file name, we don't allow reference map to be updated
        set outMapId=-1;
      END IF;
    END IF;
  ELSE
    -- server map
    select MapId, MapScale, Description into outMapId, outMapScale, outDescription from Maps where ReferenceFlag=0 and StudyId=inStudyId and Description=inDescription;
    IF no_rows_indicator THEN
      -- this is the first time we add this study map for this study
      start transaction;
      insert into Maps (StudyId, MapScale, Description, ReferenceFlag) values (inStudyId, inMapScale, inDescription, inReferenceFlag);
      select LAST_INSERT_ID() into outMapId;
      commit;
    END IF;
  END IF;
END
//
DELIMITER ;
