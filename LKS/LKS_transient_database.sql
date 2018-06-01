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
