#!/usr/bin/env python3
import sys
sys.path.append('../..')
import unittest
import numpy
import pandas
import pandas.testing
import pathogist
import pathogist.io

class IntegrationTest(unittest.TestCase):
    def setUp(self):
        true_clustering_path = 'tests/integration_tests/test_data/yersinia_true_clusterings.tsv'
        self.true_clustering = pathogist.io.open_clustering_file(true_clustering_path)
        clustering_path = 'tests/integration_tests/test_data/yersinia_final_clustering.tsv'
        self.clustering = pathogist.io.open_clustering_file(clustering_path)
       
    def test_integration(self):
        pandas.testing.assert_series_equal(self.clustering['Consensus'],self.true_clustering['Consensus'])
