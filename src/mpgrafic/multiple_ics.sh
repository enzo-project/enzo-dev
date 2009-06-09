#!/bin/sh

first=200902121
nn=10
for i in `seq 0 $nn`; do
    seed=`expr $i + $first`
    python make_ic.py $seed
done