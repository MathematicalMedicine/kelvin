-- Copyright (C) 2022 Mathematical Medicine LLC
-- 
-- This program is free software: you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by the Free
-- Software Foundation, either version 3 of the License, or (at your option)
-- any later version.
-- 
-- You should have received a copy of the GNU General Public License along
-- with this program. If not, see <https://www.gnu.org/licenses/>.
Use mysql;
Delete from user where User = '';
Update user set password = PASSWORD('RPW') where user = 'root';
Flush privileges;
