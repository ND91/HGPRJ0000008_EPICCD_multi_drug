#!/usr/bin/env bash

baseurl="http://raggr.usc.edu/Results/tmp/2018-12-4/8/"
tailurl="_CEU-FIN-GBR-IBS-TSI.csv"
projectIDs=(881352309 1142183493)

for $project in ${projectIDs[@]}; do
	echo $project
done
