import logging
import datetime as dt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster import hierarchy as hcl

logger = logging.getLogger(__name__)

def visualize(distance, name, metadata = None, columns = None, save_path = None):
    '''
    :param distance: distance matrix
    :param name: name of the sample (and name of the pdf)
    :param metadata: location of metadata file
    :param columns: list of columns in metadata file for labels on hierarchical clustering dendogram
    :param create_pdf: variable for the creation of an output pdf
    '''

    distance_histogram(distance, name, save_path)

def distance_histogram(distance, name='SAMPLE', save_path=None):
    '''
    Creates a distance histogram from the distances in the
    lower triangle of the distance matrix

    Used to find the threshold for correlation clustering
    '''

    logger.debug("Creating distance histogram")

    # Obtain distances from lower triangle
    values = distance.values
    distances = values[np.tril_indices(len(distance.index), -1)]

    plt.style.use('ggplot')
    plt.figure(0)
    n,bins,patches = plt.hist(distances, bins=500)

    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title('Distance Histogram for %s' % name)
    if save_path is not None:
        plt.savefig(save_path)
    plt.show()


def hierarchical_clustering(distance, name, metadata = None, columns = None, create_pdf = False, pdf = None):
    '''
    Runs hierarchical clustering and creates a dendrogram and MDS plot
    '''

    meta = pd.read_csv(metadata, header=0, index_col=0, sep="\t")

    # Select a subset of samples common to the distance matrix and the metadata file
    isect = [rowname for rowname in list(distance.index) if rowname in list(meta.index)]
    meta = meta.loc[isect]
    distance = distance.loc[isect, isect]

    # Run hierarchical clustering
    logger.debug("Running hierarchical clustering")
    linkage = hcl.linkage(distance, method='complete')

    # Plotting dendrogram for original labels
    logger.debug("Generating dendrograms")
    plt.figure(0)
    dendro = hcl.dendrogram(linkage, labels=distance.index)

    # Plotting dendrograms for metadata labels
    for index,col in enumerate(columns):
        plt.figure(index+1)
        dendro = hcl.dendrogram(linkage, labels=meta[col])

def plot_nawc_values(nawc_values,title):
    title = species + " " + genotype + " " + algorithm
    x_values = nawc_values.index.values.astype(numpy.float)
    y_values = nawc_values.values
    plt.title(title)
    plt.plot(x_values,y_values)
    plt.ylabel('nAWC')
    plt.xlabel('threshold')
    plt.show() 
