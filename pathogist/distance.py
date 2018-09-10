import logging
import sys
import itertools
import numpy
import pandas

logger = logging.getLogger(__name__)

def hamming_distance(calls1,calls2):
    '''
    Given two numpy 1-dimensional arrays calls1, calls2 of the same shape, returns the hamming distance.
    '''
    return numpy.count_nonzero(calls1 != calls2) 

def l1_norm(calls1,calls2):
    '''
    Given two numpy 1-dimensional arrays calls1, calls2 of the same shape, returns the L1 norm.
    '''
    return numpy.linalg.norm(calls1-calls2,ord=1)

def create_mlst_distance_matrix(calls):
    '''
    Given a dictionary of MLST calls (where sample names are keys to a vector), creates an MLST
    distance matrix represented as a Pandas Dataframe object.
    Distance: hamming distance
    '''
    samples = calls.keys()
    num_samples = len(samples)
    logger.debug("Got " + str(num_samples) + " samples...")
    distance_matrix = pandas.DataFrame(numpy.zeros(shape=(num_samples,num_samples),dtype=int),
                                       index=samples,columns=samples,dtype=int) 
    for sample1,sample2 in itertools.combinations(samples,2):
        distance_matrix[sample1][sample2] = hamming_distance(calls[sample1],calls[sample2])
        distance_matrix[sample2][sample1] = distance_matrix[sample1][sample2]
    return distance_matrix

def create_cnv_distance_matrix(calls):
    '''
    Given a dictionary of CNV calls (where sample names are keys to a vector), creates a CNV
    distance matrix represented as a Pandas Dataframe object.
    Distance: L1 norm
    '''
    samples = calls.keys()
    num_samples = len(samples)
    distance_matrix = pandas.DataFrame(numpy.zeros(shape=(num_samples,num_samples),dtype=int),
                                       index=samples,columns=samples,dtype=float) 
    for sample1,sample2 in itertools.combinations(samples,2):
        distance_matrix[sample1][sample2] = l1_norm(calls[sample1],calls[sample2])
        distance_matrix[sample2][sample1] = distance_matrix[sample1][sample2]
    return distance_matrix

def create_snp_distance_matrix(calls):
    '''
    Given a dictionary of MLST calls (where sample names are keys to a vector), creates an MLST
    distance matrix represented as a Pandas Dataframe object.
    Distance: hamming distance
    '''
    # Assume that the SNP distance matrix is created in the same was as the MLST distance matrix.
    return create_mlst_distance_matrix(calls)

def match_distance_matrices(distances):
    '''
    Modify a set of distance_matrices so that they all share the same set of samples
    '''
    samples_intersected = list(set.intersection(*[set(distances[key].columns.values)
                                               for key in distances]))
    matched_distances = {key: distances[key].loc[samples_intersected,samples_intersected]
                         for key in distances.keys()}
    return matched_distances
