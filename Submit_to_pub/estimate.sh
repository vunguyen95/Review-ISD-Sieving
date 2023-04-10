#! /bin/bash

time python3 estimate.py 2 60 2>&1 | tee > resultBIKE_HQC_2_60.txt &
time python3 estimate.py 2 120 2>&1 | tee > resultBIKE_HQC_2_120.txt &
time python3 estimate.py 5 60 2>&1 | tee > resultBIKE_HQC_5_60.txt &
time python3 estimate.py 5 120 2>&1 | tee > resultBIKE_HQC_5_120.txt 
