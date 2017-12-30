#!/bin/sh
HOST=   'host'
USER=   'user'
PASSWD= 'pass'

ftp -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASSWD
ls -1
bye
END_SCRIPT

exit 0

 
