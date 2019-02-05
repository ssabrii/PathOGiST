#!/bin/bash

set -e

./PATHOGIST run tests/integration_tests/test_data/pathogist-run_all-test.yaml 

python -m unittest tests/integration_tests/test_integration.py
# Run Genotyping Software
#./PATHOGIST run tests/integration_tests/test2_data/pathogist-run_all-test2.yaml 

#python -m unittest tests/integration_tests/test2_integration.py
./PATHOGIST run tests/integration_tests/test3_data/pathogist-run_all-test3.yaml 

python -m unittest tests/integration_tests/test3_integration.py

./PATHOGIST run tests/integration_tests/test4_data/pathogist-run_all-test4.yaml 

python -m unittest tests/integration_tests/test4_integration.py