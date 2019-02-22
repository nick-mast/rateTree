#!/bin/bash
# Script to run buildRateTree on a set of series
#Execute as "buildRateTree_batch [directory/to/rq/series']

for seriesDir in $*
do
	echo "$seriesDir/Rate.root"
	buildRateTree "$seriesDir/" "$seriesDir/Rate.root"
done

