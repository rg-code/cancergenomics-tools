#! /bin/bash

while read line ; do
wget "http://gdac.broadinstitute.org/runs/stddata__2015_02_04/samples_report/$line.html"
done < mylist.txt
