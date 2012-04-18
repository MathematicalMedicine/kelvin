CREATE TABLE Analyses (
StudyId int(11) NOT NULL,
PedigreeRegEx varchar(32) NOT NULL,
PedigreeNotRegEx varchar(32) NOT NULL DEFAULT 'XYZZY',
AnalysisId int(11) NOT NULL AUTO_INCREMENT,
InsertTime timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
PRIMARY KEY (StudyId, PedigreeRegEx, PedigreeNotRegEx),
KEY `AnalysisId` (AnalysisId)
) ENGINE=InnoDB;

CREATE TABLE RegionModels (
RegionId int NOT NULL,
ModelId int NOT NULL,
PRIMARY KEY (RegionId, ModelId),
INDEX (ModelId)
) ENGINE=InnoDB;

Alter table Regions add AnalysisId int NOT NULL default 0;
Alter table Regions drop index MinimalRegions;
Alter table Regions add unique key MinimalRegions (StudyId, AnalysisId, ChromosomeNo, RefTraitPosCM, RegionNo);
Alter table PedigreePositions add PendingLikelihoods int default 0;

Source LKS_setup_trigger_proc.sql;
