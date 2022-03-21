import numpy as np


# This is included as a seperate class to the Mapper class so it can contain more relevant information
#than a simple dictionary attribute of the mapper class
class Cluster():
    """This class run clustering on the data points contained within each interval

    Parameters
    ----------

    cover : class
        Cover class built in 'cover.py'
        Division on lens function into intervals with overlap

            """

    def __init__(self, cover):

        self.cover = cover
        self.clusters_in_intervals = {}
        self.cluster_count_in_interval = {}
        self.cluster_samples_in_interval = {}


    def run_clustering_algorithm(self, algorithm):
        """Peform clustering on the data points contained within each interval.

        Parameters
        ----------

        algorithm : clustering algorithm
            Algorithm by which to cluster the samples.
            Default is sklearn agglomerative hierarchical clustering

        """

        clusters_in_intervals = {} #each cluster in each interval set
        cluster_count_in_interval = {} #number of clusters in each set
        cluster_samples_in_interval = {}

        n = 0 #n is the number of clusters at the start
        n_i = 0 # n_i is the number of clusters generated in that interval set


        # clustering algorithms can set the minimum number of clusters
        #and will throw up an error if there are too little samples in cluster
        # n_clusters = agglo,  min_cluster_size = hdbscan. Can add more if necessary
        algorithm_parameters = algorithm.get_params()
        algorithm_minimum = 2
        options = ["n_clusters", "min_cluster_size"]
        for setting in algorithm_parameters:
            if setting in options:
                algorithm_minimum = algorithm_parameters[setting]

        #for each interval, if there is more than two data points in that interval
        #Fit a clustering algorithm to those datapoints
        for i in range(0, self.cover.intervals):
            #access the unique samples within each interval
            samples = self.cover.samples_in_interval[i]
            points = self.cover.data_in_interval[i]

            #THIS SETS THE MIMINMUM NUMBER OF SAMPLES IN A NODE - we lose intervals if we do not incorporate nodes
            if len(samples) > algorithm_minimum:
                cl = algorithm.fit(points)

                #j is each cluster in interval set i
                j_cl = [points[np.where(cl.labels_ == label)] for label in set(cl.labels_)]
                clusters_in_intervals[i] = j_cl

                # n_i is the number of clusters generated in that interval set
                n_i = len(j_cl)
                cluster_count_in_interval[i] = n_i

                #c_i is the index for each cluster
                #obtain the samples from the cluster and not the index because
                #outlier samples can be dropped from the final graph model
                c_i = [samples[np.where(cl.labels_ == label)] for label in set(cl.labels_)]
                cluster_samples_in_interval[i] = c_i

            n =  n + n_i


        self.cluster_samples_in_interval = cluster_samples_in_interval
        self.cluster_count_in_interval = cluster_count_in_interval
        self.clusters_in_intervals = clusters_in_intervals
