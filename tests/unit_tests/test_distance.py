import numpy
import pandas as pd
import pandas.testing as pt
import sys
sys.path.append('../..')
import pathogist
import pathogist.io
import pathogist.distance
from unittest import TestCase

class DistanceTest(TestCase):

    def setUp(self):
        self.mlst_calls = pathogist.io.read_mlst_calls('test_data/distance/mlst_calls_path.txt')
        self.cnv_calls = pathogist.io.read_cnv_calls('test_data/distance/cnv_calls_paths.txt')

    def test_hamming_distance(self):
        vector1 = numpy.array([0,0,0,0,0,0],dtype=int)
        vector2 = numpy.array([1,0,0,1,1,0],dtype=int)
        true_distance = 3
        self.assertEqual(pathogist.distance.hamming_distance(vector1,vector2),true_distance)

    def test_l1_norm(self):
        vector1 = numpy.array([10,20,30,40,50],dtype=int)
        vector2 = numpy.array([0,10,20,30,40],dtype=int)
        true_distance = 50
        self.assertEqual(pathogist.distance.l1_norm(vector1,vector2),true_distance)

    def test_mlst_distance_matrix(self):
        true_matrix = pathogist.io.open_distance_file('test_data/distance/mlst_distance_matrix.tsv') 
        pt.assert_frame_equal(pathogist.distance.create_mlst_distance_matrix(self.mlst_calls),
                              true_matrix)

    def test_snp_distance_matrix(self):
        '''
        The SNP distance matrix creation function should work the same way as the MLST one, for now.
        '''
        true_matrix = pathogist.io.open_distance_file('test_data/distance/mlst_distance_matrix.tsv') 
        pt.assert_frame_equal(pathogist.distance.create_mlst_distance_matrix(self.mlst_calls),
                              true_matrix)

    def test_cnv_distance_matrix(self):
        true_matrix = pathogist.io.open_distance_file('test_data/distance/cnv_distance_matrix.tsv')
        pt.assert_frame_equal(true_matrix,pathogist.distance.create_cnv_distance_matrix(self.cnv_calls)) 
