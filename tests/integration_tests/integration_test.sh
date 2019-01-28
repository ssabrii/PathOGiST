#!/bin/bash

set -e

./PATHOGIST run tests/integration_tests/test_data/pathogist-run_all-test.yaml 

python -m unittest tests/integration_tests/test_integration.py
