#!/bin/bash

set -e

./PATHOGIST run tests/integration_tests/test_data/pathogist-run_all-test.yaml 

python -m unittest tests/integration_tests/test_integration.py

./PATHOGIST run tests/integration_tests/test2_data/pathogist-run_all-test2.yaml 

python -m unittest tests/integration_tests/test2_integration.py
