#!/bin/bash

for rep_no in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ; do
	mkdir -p replica_$rep_no
	cp basefiles/* replica_$rep_no
		
done
