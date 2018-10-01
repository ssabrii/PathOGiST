import numpy
import pandas as pd
import pandas.testing as pt
import sys
import unittest 
sys.path.append('../..')
import pathogist
import pathogist.io
import pathogist.cluster

class ClusterTest(unittest.TestCase):

    def setUp(self):
        mlst_dist_path = 'tests/unit_tests/test_data/cluster/yersinia_mlst_dist.tsv'
        snp_dist_path = 'tests/unit_tests/test_data/cluster/yersinia_snp_dist.tsv'
        kwip_dist_path = 'tests/unit_tests/test_data/cluster/yersinia_kwip_dist.tsv'
        self.snp_dist = pathogist.io.open_distance_file(snp_dist_path)
        self.mlst_dist = pathogist.io.open_distance_file(mlst_dist_path)
        self.kwip_dist = pathogist.io.open_distance_file(kwip_dist_path)

        mlst_clust_path = 'tests/unit_tests/test_data/cluster/yersinia_mlst_clustering.tsv'
        snp_clust_path = 'tests/unit_tests/test_data/cluster/yersinia_snp_clustering.tsv'
        kwip_clust_path = 'tests/unit_tests/test_data/cluster/yersinia_kwip_clustering.tsv'
        self.mlst_clust = pathogist.io.open_clustering_file(mlst_clust_path)
        self.snp_clust = pathogist.io.open_clustering_file(snp_clust_path)
        self.kwip_clust = pathogist.io.open_clustering_file(kwip_clust_path)

        self.clusterings = {'SNP': self.snp_clust, 'MLST': self.mlst_clust, 'KWIP': self.kwip_clust}
        self.distances = {'SNP': self.snp_dist, 'MLST': self.mlst_dist, 'KWIP': self.kwip_dist}
        self.fine_clusterings = ['SNP']

    def test_correlation(self):
        true_clustering_path = 'tests/unit_tests/test_data/cluster/yersinia_true_clustering.tsv'
        true_clustering = pathogist.io.open_clustering_file(true_clustering_path)
        clustering = pathogist.cluster.correlation(self.mlst_dist,500)
        samples = list(clustering.index.values)
        self.assertEqual(pathogist.cluster.adjusted_rand_index(clustering.loc[samples],
                                                               true_clustering.loc[samples]),1)

    def test_consensus(self):
        clustering = pathogist.cluster.consensus(self.distances,self.clusterings,self.fine_clusterings)
        true_clustering_path = 'tests/unit_tests/test_data/cluster/yersinia_consensus_clustering.tsv'
        true_clustering = pathogist.io.open_clustering_file(true_clustering_path)
        samples = list(clustering.index.values)
        pt.assert_series_equal(clustering.loc[samples,'Consensus'],
                              true_clustering.loc[samples,'Consensus'])
