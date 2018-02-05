Create database LKS;
use mysql;
Delete from user where User = '';
Grant all on *.* to 'adminLKS'@'localhost' identified by 'barfewww';
Grant all on *.* to 'adminLKS'@'%' identified by 'barfewww';
flush privileges;
