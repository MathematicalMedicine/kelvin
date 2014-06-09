Use mysql;
Delete from user where User = '';
Update user set password = PASSWORD('RPW') where user = 'root';
Flush privileges;
