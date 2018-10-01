#!/bin/bash

set -e

../../pathogist.py all test_data/pathogist-run_all-test.yaml

python -m unittest test_integration.py
