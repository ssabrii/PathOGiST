#!/bin/bash

set -e

./PATHOGIST all tests/integration_tests/test_data/pathogist-run_all-test.yaml 

python -m unittest tests/integration_tests/test_integration.py
