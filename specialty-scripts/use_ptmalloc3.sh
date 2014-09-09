#!/bin/bash

case $HOSTNAME in
    Levi-Montalcini* )
	echo "yes"
#        WERROR=-Werror
        ;;
    * )
	echo "no"
        WERROR=
	;;
esac
