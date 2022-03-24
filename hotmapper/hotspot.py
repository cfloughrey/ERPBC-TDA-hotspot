# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 10:57:31 2020

@author: ciara
"""

import numpy as np
import networkx as nx
import pandas as pd
import hotmapper.utils as hmu

#plotting
import matplotlib as mpl
import matplotlib.pyplot as plt

#dendogram
from scipy.spatial import distance
from scipy.cluster import hierarchy

plt.rcParams.update({'font.size': 22})




#----------------------hotspot class--------------------------------#
class HotSpot:
    """ The hotspot class searches a collection of interconnected nodes
    obtained from the Mapper grah, identifying any clusters of nodes that
    present anomalous attribute filter_values

    Parameters
    ----------

    graph : networkx graph
        Networkx graph obtained from the Mapper class

    attribute : int or float list
        Attribute function of predefined values for each node to indicate the colouring of the nodes
        """


    def __init__(self, mapper, subgraph, subgraph_number, text = True):

        #each subgraph needs to be assigned new labels as an new graph, previously contain
        #labels for original mapper. labels retained in nx attributes'_node'
        self.mapper = mapper
        self.old_labels = subgraph
        self.graph = nx.convert_node_labels_to_integers(subgraph, label_attribute='old_label')
        self.subgraph_number = subgraph_number

        self.edge_cutoff = None
        self.Z = [] #linkage matrix
        self.edge_weight_sorted = {}


    def identify_cutoffpoint(self):
        """ Define edge weights according to the difference in attribute value between nodes.
            Sort the edge weights, then identify the cut-off point to build the clusters - between the
            edge weights with the largest difference in values"""

        #assigns edge weights for every pair of edges in the original graph
        edge_w_org = []
        for e in self.old_labels.edges():
            #find the absolute difference in nodes
            abs_diff = abs(self.mapper.attribute[e[0]] - self.mapper.attribute[e[1]])
            self.old_labels[e[0]][e[1]]['weight'] = abs_diff #assign as edge weight
            edge_w_org.append(abs_diff)

        #create a dictionary containing the correct weights for new labels
        edge_w = {(e[0],e[1]): edge_w_org[i] for i,e in  enumerate(self.graph.edges())}

        #sort labels in order of ascending weight
        ew_sort = {k: v for k, v in sorted(edge_w.items(), key=lambda item: item[1])}

        #empty matrix of node length x needed for linkage tree
        a = np.ones(shape=(len(self.graph.nodes),len(self.graph.nodes)))

        #construct matrix of nodes and edges
        for i,j in ew_sort.items():
            u = i[0]
            v = i[1]
            a[u,v] = j
            a[v,u] = j
            a[u,u] = 0
            a[v,v] = 0

        #construct a single linkage matrix of connected nodes
        dists = distance.squareform(a)
        Z = hierarchy.linkage(dists)

        #retain the threshold for edge cutoff and get the edge weight values from the dict
        #find the differences between all the values. We ignore the very last value
        l = list(ew_sort.values())
        diff = [j-i for i, j in zip(l[:-2], l[1:])]

        #find the maximum difference and the index
        try:
            m = max(diff)
            max_index = diff.index(m)
            cut_index = max_index
        except ValueError:
            m = 0
            cut_index = 0
        #index at the cutoff value + half way between the difference
        cutoff = l[cut_index]+ m*0.5

        self.edge_cutoff = cutoff
        self.Z = Z #linkage matrix
        self.edge_weight_sorted = ew_sort



    def cluster_identification(self):
        """Function seperates the graph into clusters according to the edge
        weight cutoff point. Any edges above this threshold are removed, leaving
        a graph seperated into groups of similar node values"""

        #create a new graph object for each cluster
        H = nx.convert_node_labels_to_integers(self.graph)

        #obtain the weight values as list
        w_s = [v for i,v in self.edge_weight_sorted.items()]

        #if edge weight cutoff is supplied
        w_s_cutoff = [i for i in w_s if i <= self.edge_cutoff]
        r = len(self.edge_weight_sorted) - len(w_s_cutoff)
        remove_edges = list(self.edge_weight_sorted.items())[-r:]

        #remove these values from the graph
        for k,v in remove_edges:
            H.remove_edge(k[0],k[1])

        #map between the old and new labels
        mapping = dict(zip(self.graph.nodes, self.old_labels.nodes))
        graph_clusters = nx.relabel_nodes(H, mapping, copy = False)

        communities = hmu.obtain_community_partition(graph_clusters)
        cluster_nodes = [list(cluster.nodes) for cluster in list(communities.values())]

        return cluster_nodes




    def cluster_classification(self, clusters, extreme, attribute_threshold, min_sample_size):
        """Each identified cluster within the subgraph is classified as a hotspot or not
        according to the minimum sample size and average attribute value compared to the
        neighbourhood nodes
            """

        ## Set up the hotspot dictionary containing information for each subgraph ##
        hotspot = {}

        #find all the nodes in this graph community
        hotspot["component"] = set(hmu.flatten(clusters))
        hotspot["clusters"] = clusters
        hotspot["mean attribute"] = hmu.community_filterval(hotspot["clusters"], self.mapper.attribute)
        #Size: Nodes have overlapping patients. If you sum each node size the total will be greater than the cohort
        #the solution is to find the sample size for each cluster of nodes
        hotspot["size"] = hmu.community_size_samples(hotspot["clusters"], self.mapper.samples_in_node)
        #define boolean values for whether a hotspot is a local and /or global hotspo
        hotspot["hotspot class"] = [True]*len(hotspot["clusters"])
        # Neighbourhood: All remaining nodes from hotspot
        hotspot["neighbourhood nodes"] = hmu.build_neighbourhood(hotspot)
        hotspot["neighbourhood size"] = hmu.community_size_samples(hotspot["neighbourhood nodes"], self.mapper.samples_in_node)
        hotspot["neighbourhood attribute"] = hmu.community_filterval(hotspot["neighbourhood nodes"], self.mapper.attribute)



        ## CHECK 1 - Size of samples in the cluster is sufficiently large ##
        #if the size of the cluster is less than a predefined threshold
        #classify it as a non-hotspot and remove it from C
        for i,c in enumerate(hotspot["size"]):
            if c < min_sample_size:
                hotspot["hotspot class"][i]=False


        ## CHECK 2 - Size of samples in the cluster is smaller than neighbourhood ##
        #calculate the threshold of the clusters
        #m1 is the threshold based on median absolute deviation for the community sizes
        m1 = hmu.mad(hotspot["size"])
        #if the S(Nc) - S(C) < MAD threshold then classify as non hotspot
        for i, cluster_size in enumerate(hotspot["size"]):
            neighbour_size = hotspot["neighbourhood size"][i]
            if (neighbour_size - cluster_size) < m1:
                hotspot["hotspot class"][i] = False


        ## CHECK 3 - attribute difference between hotspot and neighbourhood is sufficient ##
        #If abs(F(C) - F(Nc) < threshold then consider hotspot as False
        for i, cluster_att in enumerate(hotspot["mean attribute"]):
            #find the neighbourhood attribute for that cluster
            neighbour_att = hotspot["neighbourhood attribute"][i]

            #allow the user to specify which range of values to search for
            #e.g. user specifies "lower". If the difference between attribute and cluster
            #is above the threshold BUT the cluster is higher than the neighbour, set lower check to True
            #and reject hotspot
            low_check = (neighbour_att < cluster_att)
            high_check = (neighbour_att > cluster_att)
            att_check = abs(neighbour_att - cluster_att) < attribute_threshold

            extreme_options = {"lower": any([att_check, low_check]),
                                "higher": any([att_check, high_check]),
                                "either": att_check}

            if extreme_options[extreme] == True:
                hotspot["hotspot class"][i] = False

        #hotspot dataframe is returned regardless if a hotspot is found to be true or false
        hotspot["hotspot nodes"] = [hotspot["clusters"][i] for i, x in enumerate(hotspot["hotspot class"]) if x]
        hotspot["hotspot samples"] = hmu.community_samples(hotspot["hotspot nodes"], self.mapper.samples_in_node)
        return hotspot




    def plot_dendrogram(self,  labels = None):
        """Plot the dendogram illustrating at what edge weight value do nodes become
        connected. Nodes are coloured by attribute value"""

        #obtain the hex colours of nodes to assign colours to labels
        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=min(self.mapper.attribute), vmax=max(self.mapper.attribute))
        col_hex = [mpl.colors.to_hex(c, keep_alpha=False) for c in cmap(norm((self.mapper.attribute)))]

        #plot figure
        plt.clf()
        plt.figure()
        plt.rc('font', size=22)          # controls default text sizes
        plt.rc('axes', titlesize=22)     # fontsize of the axes title
        plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=16)
        plt.rc('legend', fontsize=16) # fontsize of the tick labels
        plt.rc('figure', titlesize=20)  # fontsize of the figure title

        plt.title(F"Connectivity of nodes: subgraph {self.subgraph_number}")
        dn = hierarchy.dendrogram(self.Z, color_threshold= self.edge_cutoff, above_threshold_color='grey', labels = labels)
        ax = plt.gca()
        ax.tick_params(axis='x', labelsize=16)
        xlbls = ax.get_xmajorticklabels()
        for lbl in xlbls:
            lbl.set_color(col_hex[int(lbl.get_text())])
        plt.show()

        return dn









#--------------- RUN HOTSPOT DETECTION FOR EACH COMPONENT ----------------------
class Subgraphs:
    """This class contains all the connected component subgraphs of the original
     mapper graph as seperate entities. Hotspot detection will be run on each component

    Parameters
    ----------

    graph : networkx graph
        Networkx graph obtained from the Mapper class

    attribute : int or float list
        Attribute function of predefined values for each node to indicate the colouring of the nodes
        """

    def __init__(self, mapper, text = True):

        #We seperate the mapper class into the interconnected subgraphs (S)
        self.mapper = mapper
        self.subgraphs = hmu.obtain_community_partition(mapper.graph)
        self.text = text

    def run_hotspot_search(self, attribute_threshold, min_sample_size = 30, extreme = "either", dendrogram = False):
        """run a search for hotspots of nodes on each connected subgraph.

        Parameters
        ----------

        attribute_threshold : int, float
            Value specifying the threshold at which to accept a hotspot. If the difference in
            attribute value between a group of nodes is larger than this (absolute) threshold
            in comparison to their neighbours, then accept as a hotspot.

        extreme : ["lower","higher", "none"], default: ```either```
            Whether to specify which extremity of the attribute function to search for hotspots
            Default is either high or low attribute hotspot

        min_sample_size : int, default: ```30```
            The minimum number of samples within a group of nodes to be confirmed as a hotspot.
            Default is 30

        dendogram : boolean, default: ```False```
            If set to true, plots a the dendrogram of connections between nodes seperated
            by attribute values


        Returns
        ----------
        hotspots : dictionary
            Returns a dictionary of hotspots found within each subgraph
            """

        if self.text == True:
            print("Searching for hotspots within the subgraphs...")

        #Run hotspot search within each subgraph if the subgraph has more than 2 nodes
        hotspot_collection = {}
        for i in self.subgraphs:
            graph = self.subgraphs[i]
            node_no = len(list(graph.nodes))
            subgraph_number = i

            if node_no > 2:
                #initialise hotspot class
                hotspot = HotSpot(self.mapper, graph, subgraph_number)

                # Step 1 - identify the hotspots in the graph
                hotspot.identify_cutoffpoint()
                clusters = hotspot.cluster_identification()

                #Step 2 - Classify clusters as hotspot or not
                hotspot_clusters = hotspot.cluster_classification(clusters,
                                                                  extreme,
                                                                  attribute_threshold,
                                                                  min_sample_size = int(min_sample_size))
                hotspot_collection[i] = hotspot_clusters

                #plot dendogram if specified
                if dendrogram == True:
                    hotspot.plot_dendrogram(labels = list(graph.nodes))

            else:
                pass

        #return the hotspot nodes in the graph
        if self.text == True:
            for i in hotspot_collection:
                check = hotspot_collection[i]['hotspot class']
                if any(check):
                    print(F"Hotspots found in subgraph {i}: nodes {hotspot_collection[i]['hotspot nodes']}")
                else:
                    print(F"No hotspots present in subgraph {i}")

        hotspot_df = pd.DataFrame(hotspot_collection)
        return hotspot_df
