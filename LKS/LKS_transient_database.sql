Create database DB;
Use mysql;
Grant all on DB.* to 'USER'@'localhost';
Grant all on DB.* to 'USER'@'%';
Update user set password = PASSWORD('UPW'),SUPER_PRIV='Y',FILE_PRIV='Y',PROCESS_PRIV='Y' where user='USER';
Flush privileges;
