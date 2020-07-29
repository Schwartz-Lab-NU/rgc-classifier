#!/bin/bash

for file in `ls *params.in`;
do ./train "${file:0:4}" "10">> "${file:0:4}.log" 2>&1;
done
