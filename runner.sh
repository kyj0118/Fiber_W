#!/bin/bash

mm=2
#mm=1
#mm=0.5

for (( en=100; en<2100; en+=100 ))
#for (( en=100; en<200; en+=100 ))
do
    bsub -q s ./exampleB5 run_gamma${en}MeV.mac ../root/gamma${en}MeV_${mm}mm 1
done
