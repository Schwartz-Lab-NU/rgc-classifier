#!/bin/bash

for file in `ls *params.in`;
do ./ecoc "${file:0:4}" >> "${file:0:4}.log" 2>&1;
done
