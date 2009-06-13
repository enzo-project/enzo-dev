#!/bin/sh

for i in *.f90; do
    base=${i%%.*}
    mv $i /tmp/${base}.F90
    mv /tmp/${base}.F90 ${base}.F90
done