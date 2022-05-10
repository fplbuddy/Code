#!/bin/bash

val=0.0001
strin=E${val}
destin=/travail/skannabiran/Vorder/Size1024/Run_0.0000000625_0.5/

make dist clean mhd2DB
mv mhd2DB mhd2DB_run
cp {status.txt,parameter.txt} ${destin}${strin}
cp mhd2DB_run ${destin}${strin}
clear

