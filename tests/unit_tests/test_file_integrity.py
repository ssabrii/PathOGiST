import numpy
import pandas as pd
import pandas.testing as pt
import sys
import unittest 
sys.path.append('../..')
import pathogist
import pathogist.io
import pathogist.cluster
import yaml

class FileIntegrityTest(unittest.TestCase):

    def setUp(self):
        config_path = 'tests/unit_tests/test_data/file_integrity/config.yaml'


    def test_config(self):
        with open(self.config_path,'r') as config_stream:
            try:
                config = yaml.load(config_stream) 
            except yaml.YAMLError:
                print(yaml.YAMLError)
                sys.exit(1)
        assert pathogist.assert_config(config) == 0

    '''
    def test_fastq_input(self):
        clustering = pathogist.cluster.consensus(self.distances,self.clusterings,self.fine_clusterings)
        true_clustering_path = 'tests/unit_tests/test_data/cluster/yersinia_consensus_clustering.tsv'
        true_clustering = pathogist.io.open_clustering_file(true_clustering_path)
        samples = list(clustering.index.values)
        pt.assert_series_equal(clustering.loc[samples,'Consensus'],
                              true_clustering.loc[samples,'Consensus'])

    '''