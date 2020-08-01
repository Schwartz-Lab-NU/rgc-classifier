#!/bin/bash

for file in `ls *params.in`
do
	echo "Starting paramset ${file:0:4}";
	./ecoc "${file:0:4}" >> "${file:0:4}.log" 2>&1;
	echo -e "\tTraining...";
	for i in {1..3}
	do
		./train "${file:0:4}" "15">> "${file:0:4}.log" 2>&1;
	done
	echo -e "\tTransforming...";
	for i in {1..9}
	do
		./transform "${file:0:4}">> "${file:0:4}.log" 2>&1;
	done
	echo -e "\tEvaluating...";
	
	for i in {1..3}
	do
		for j in {1..3}
		do
			for k in {1..3}
			do
				./test "${file:0:4}" "$i" "$j" "$k">> "${file:0:4}.log" 2>&1;
			done
		done
	done

	echo "Done paramset ${file:0:4}";
done
