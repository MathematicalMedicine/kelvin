# Grab STUDY directive for database parameters
study=$(grep -i ^Study client.conf)
set -- $study
Status=$(mysql --host $4 --user $6 --password=$7 $5 --batch --skip-column-names --execute="call CleanStudy('$2');")
