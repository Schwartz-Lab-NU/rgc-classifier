#!/bin/bash

for file in `ls *params.in`
do
	echo "Starting paramset ${file:0:4}";
	./train "${file:0:4}" "15" >> "${file:0:4}.log" 2>&1;
	echo "Done paramset ${file:0:4}";
done
