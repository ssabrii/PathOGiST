#!/usr/bin/env python3
import unittest
import numpy
import pandas
import pandas.testing
sys.path.append('../..')
import pathogist
import pathogist.io

def IntegrationTest(unittest.TestCase):
    def setUp(self):
        true_clustering_path = 'tests/integration_tests/test_data/yersinia_true_clustering.tsv'
        self.true_clustering = pathogist.io.open_clustering_file(true_clustering_path)
        clustering_path = 'tests/integration_tests/test_data/yersinia_final_clustering.tsv'
        self.clustering = pathogist.io.open_clustering_file(clustering_path)
       
    def test_integration(self):
        pandas.testing.assert_frame_equal(self.clustering,self.true_clustering)
