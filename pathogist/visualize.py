import logging
import pandas 
import numpy 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster import hierarchy as hcl
import sklearn
import sklearn.manifold
import networkx
import random
import itertools

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
    distances = values[numpy.tril_indices(len(distance.index), -1)]

    plt.style.use('ggplot')
    plt.figure(0)
    n,bins,patches = plt.hist(distances, bins=500)

    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.title('Distance Histogram for %s' % name)
    if save_path is not None:
        plt.savefig(save_path)
    plt.show()

def visualize_clusterings(summary_clustering):
    '''
    Visualize the consensus clustering and correlation clusterings.
    @param summary_clustering: A Pandas Dataframe where the index is the sample names, columns are
                               clusterings, and entries are cluster assignments.
    '''
    samples = summary_clusterings.index.values
    sample_positions = {sample: (random.randrange(0,100),random.randrange(0,100)) for sample in samples}
    clusterings_names = summary_clusterings.columns.values
    
    for clustering_name in clustering_names:
        plt.figure()
        graph = networkx.Graph()
        for sample in samples:
            graph.add_node(sample,Position=sample_positions[sample]) 
        clustering = summary_clusterings[clustering_name]
        num_clusters = numpy.amax(clustering.values)
        for cluster in range(1,num_clusters+1):
            cluster_samples = clustering.index[clustering[clustering_name] == cluster].tolist()
            for sample1,sample2 in itertools.combinations(cluster_samples,2)
                graph.add_edge(sample1,sample2)
        networkx.draw(graph)
        plt.show()
        
def hierarchical_clustering(distance, name, metadata = None, columns = None, create_pdf = False, pdf = None):
    '''
    Runs hierarchical clustering and creates a dendrogram and MDS plot
    '''

    meta = pandas.read_csv(metadata, header=0, index_col=0, sep="\t")

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
