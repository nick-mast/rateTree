#!/bin/bash
# Script to run buildRateTree on a set of series

for seriesDir in $*
do
	echo "$seriesDir/Rate.root"
	buildRateTree "$seriesDir/" "$seriesDir/Rate.root"
done

