# Copyright (C) 2010, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

# Grab STUDY directive for database parameters
study=$(grep -i ^Study client.conf)
set -- $study
Status=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call CleanStudy('$2');")
