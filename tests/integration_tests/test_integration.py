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
        self.true_clustering.sort_index(axis=1)
        clustering_path = 'tests/integration_tests/test_data/yersinia_final_clustering.tsv'
        self.clustering = pathogist.io.open_clustering_file(clustering_path)
        self.clustering.sort_index(axis=1)
       
    def test_integration(self):
        pandas.testing.assert_frame_equal(self.clustering.sort_index(axis=1),self.true_clustering.sort_index(axis=1))
