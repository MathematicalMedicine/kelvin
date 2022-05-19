-- Copyright (C) 2022 Mathematical Medicine LLC
-- 
-- This program is free software: you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the Free
-- Software Foundation, either version 3 of the License, or (at your option)
-- any later version.
-- 
-- You should have received a copy of the GNU General Public License along
-- with this program. If not, see <https://www.gnu.org/licenses/>.

CREATE DATABASE _DB_;
-- The 'WITH mysql_native_password' bit is because we still have pre-8.0
-- clients that have to be able to connect to our 8.0 server
-- also, 8.0 wants us to explicitly create a user before we grant it privs
CREATE USER '_USER_'@'%' IDENTIFIED WITH mysql_native_password BY '_UPW_';
GRANT ALL ON _DB_.* TO '_USER_'@'%';
-- The privs granted below seem to have been for debugging purposes only;
-- either they're no longer necessary (FILE) or we usually just root anyways 
-- (PROCESS, SUPER).
--GRANT SUPER,FILE,PROCESS ON *.* TO '_USER_'@'%';
FLUSH PRIVILEGES;
