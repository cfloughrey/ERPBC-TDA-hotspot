
import numpy as np
import sklearn.cluster as sklc
from sklearn import manifold, decomposition
import networkx as nx

#plotting
import matplotlib as mpl
import matplotlib.pyplot as plt

#supporting python scripts
from hotmapper.DSGA_transformation import DSGA
import hotmapper.feature_combination_lens as hm_lens
import hotmapper.utils as hm_utils
from hotmapper.covering import Cover
from hotmapper.clustering import Cluster

#single linkage dendogram
from scipy.spatial import distance

plt.rcParams.update({'font.size': 22})





class Mapper():
    """ The Mapper class builds a network graph from data. The workflow allows you to:
            1. Transform the dataset
            2. Construct a lens function to project the data to a lower dimension
            4. Build a cover of the lens from overlapping intervals
            3. Perform clustering within each interval
            5. Construct a graph from the clustering
            6. Visualise the graph using networkx

        Parameters
        ----------

        X : numpy array
            The data to build the graph on
            """


    def __init__(self, X, text = True):

        self.data = X
        self.text = text

        self.lens = None
        self.cover = None
        self.clustering = None
        self.random_value = []
        self.random_index = []
        self.node_dict = {}
        self.samples_in_node = {}
        self.samplecount_node = {}
        self.graph = None

        print("Initializing Mapper class...")



    def lens_function(self, selection = None, predefined_lens = None, g = 0.5, r = None, nr = None):
        """Project the data according to the lens function.

        Parameters
        ----------

        selection : str, int, default: ``mean``
            Can be one of several numpy statistical methods ["sum", "mean", "median", "max", "min", "std"] or euclidean distance ["l2norm"]
            Or sklearn dimensionality reduction model. CUrrent options are ["tsne", "pca"]
            Can be random weighted combination of features - To incorporate all features specificy one of ["linear", "non-linear"]. For a subset of features specify
                one of ["linear_subset", "non_linear_subset"]. If subset of features define parameter 'g'

        g : float
            Specify the fraction of features to subset the random weighted combination of feature lens function (e.g. 0.5 = 50% of features used to build lens)

        Returns
        -------
        lens : Numpy Array
            Transformed data

        """
        if self.text == True:
            print(F"Projecting data across {selection} lens...")

        #predetermined lens
        try:
            if predefined_lens is not None:
                self.lens = predefined_lens
                return predefined_lens
        except ValueError:
            raise ValueError("If input is own lens options, size of lens must equal number of samples")

        #standard numpy operations
        standard_lens = {"sum": np.sum,
                        "mean": np.mean,
                        "median": np.median,
                        "max": np.max,
                        "min": np.min,
                        "std": np.std}

        # dimensionality reduction models
        dimred_lens = {"tsne": manifold.TSNE(n_components=2, random_state=42).fit(self.data),
                        "pca": decomposition.PCA(n_components=2, random_state=42).fit(self.data)}

        #random weighted combinations of features
        weighted_features = {"linear": hm_lens.linear,
                             "non_linear": hm_lens.non_linear,
                             "linear_subset": hm_lens.linear_subset,
                             "non_linear_subset": hm_lens.non_linear_subset}

        #if lens is list of numbers, return the dimensions specified
        if type(selection) == list:
            lens = self.data[:,selection]
            self.lens = lens
            return lens

        #try the input against the lens options provided
        #euclidean norm
        elif selection == "l2norm":
            lens = hm_utils.norm(2, 1, self.data)
            self.lens = lens
            return lens

        elif selection in standard_lens:
            lens = standard_lens[selection](self.data, axis=1)
            self.lens = lens
            return lens

        elif selection in dimred_lens:
            dr = dimred_lens[selection]
            lens = dr.fit_transform(self.data)
            self.lens = lens
            return lens

        elif "subset" in selection:
            lens, r, nr = weighted_features[selection](self.data, g,r,nr)
            self.random_value = r
            self.random_index = nr
            self.lens = lens
            return lens

        else:
            try:
                lens, r = weighted_features[selection](self.data, r)
                self.random_value = r
                self.lens = lens
                return lens

            except ValueError:
                raise ValueError("Lens must be valid option. See documentation for more information")



    def covering(self, intervals = 10, overlap = 0.2):
        """Build a cover of the lens function. Currrently this implementatation only accepts 1-D lens to build the cover.
        The cover deconstructs the lens into overlapping intervals.

        Parameters
        ----------

        intervals : int, default: ``10``
            Number of sections to divide the lens function into

        overlap : float, default: ``0.2``
            Fraction of overlap between the intervals of the lens function

        Returns
        -------

        cover : dictionary
            Keys contain the split of lens into interval range
            Values contain the original data points contained in each interval

            """
        if self.text == True:
            print(F"Building cover over lens of {intervals} intervals with {overlap * 100}% overlap")

        #The mapper class function accesses the Cover() class in a seperate script
        #this allows me to contain all the attributes for the cover together
        cover =  Cover(self.data, self.lens, intervals, overlap)
        cover.build_cover()

        self.cover = cover



    def cluster_data(self, algorithm = sklc.AgglomerativeClustering()):
        """Peform clustering on the data points contained within each interval.

        Parameters
        ----------

        algorithm : clustering algorithm, default: ``sklc.AgglomerativeClustering()``
            Algorithm by which to cluster the samples.
            Default is sklearn agglomerative hierarchical clustering

            """
        if self.text == True:
            print(F"Run clustering algorithm '{algorithm}' within across the cover...")

        #The mapper class function accesses the Cover() class in a seperate script
        #this allows me to contain all the attributes for the cover together
        clusterclass =  Cluster(self.cover)
        clusterclass.run_clustering_algorithm(algorithm)

        self.clustering = clusterclass



    def build_graph(self, attribute):
        """Build a graph on the dataset where a node is constructed for each
        cluster and edges formed between clusters with overlapping samples


        Parameters
        ----------
        attribute : int or float list
            Attribute function of predefined values for each node to indicate the colouring of the nodes

        Returns
        -------

        G : networkx graph
            Returns an undirected network graph. Nodes represent clusters of samples and edges represent overlapping clusters

            """
        if self.text == True:
            print("Building graph...")

        #networkx graph class
        G = nx.Graph()

        node = 0
        node_dict = {}
        samples_in_node = {}
        samplecount_node = {}

        #construct a graph with a vertex for each cluster
        #access the samples in each cluster in each interval to build a node
        for i in self.clustering.cluster_count_in_interval:
            #for each cluster in that interval
            for j in range(0, self.clustering.cluster_count_in_interval[i]):
                #define the node index values simultaneously
                samples_in_node[node] = self.clustering.cluster_samples_in_interval[i][j]
                samplecount_node[node] = len(np.unique(self.clustering.cluster_samples_in_interval[i][j]))
                #add the node of samples to the graph
                G.add_node(node)
                #build the node dictionary to understand which intervals and clusters the nodes correspond to
                node_dict[(i,j)] = node
                node = node + 1


        #to build an edge check for overlapping samples between neighbouring nodes
        #define the key list to iterate through dict and check neighbouring intervals
        keyList=sorted(self.clustering.cluster_samples_in_interval.keys())
        for i, l in enumerate(keyList[:-1]):
            l_neighbour = keyList[i+1]
            cluster_points = self.clustering.cluster_samples_in_interval[l]
            cluster_points_neighbour = self.clustering.cluster_samples_in_interval[l_neighbour]

            for cluster_label, main_cluster in enumerate(cluster_points):
                for cluster_label_neighbour, neighbour_cluster in enumerate(cluster_points_neighbour):
                    #if any points in the main cluster being check are in the neighbouring cluster
                    if any(point in main_cluster for point in neighbour_cluster):
                        G.add_edge(node_dict[(l,cluster_label)], node_dict[(l_neighbour,cluster_label_neighbour)])

        #assign predefined list for node colour if specified
        #else colour by lens
        try:
            self.attribute = hm_utils.colour_by_y(samples_in_node, attribute)
        except ValueError:
            print("Please specify the colouring of the graph")


        self.node_dict = node_dict
        self.samples_in_node = samples_in_node
        self.samplecount_node = samplecount_node
        self.graph = G

        return G




    def visualise(self, style = 1, size = 10, labels = False, col_legend_title = "Legend"):
        """Visualise the networkx graph.

        Parameters
        ----------

        style : [1],[2], default: ``1``
            Selects either fruchterman_reingold_layout from networkx (1) or kamada_kawai_layout (2) to structure the graph

        size : int, default: ``10``
            Size of node legends specifying the number of samples per nodes

        labels = boolean, default: ``False``
            Specifies whether to label nodes with number

        col_legend_title = str, default: ```Legend```
            Labels the attribute legend

            """

        if self.text == True:
            print("Visualising graph...")

        network_styles = {1: nx.fruchterman_reingold_layout(self.graph, seed=300),
                          2: nx.kamada_kawai_layout(self.graph)}

        #Plot figure with legend, specifying style
        cmap = mpl.cm.viridis
        fig = plt.figure(figsize=(12, 12), constrained_layout=True)
        pos =  network_styles[style]
        norm = mpl.colors.Normalize(vmin=min(self.attribute), vmax=max(self.attribute))

        #specify the number of samples in each node according to size attribute
        nsize = hm_utils.node_size(self.samples_in_node)
        if size == 0:
            nodes = nx.draw_networkx_nodes(self.graph,
                                          pos= pos,
                                          node_color=cmap(norm((self.attribute))),
                                          alpha=1)
        else:
            nodes = nx.draw_networkx_nodes(self.graph,
                                          pos= pos,
                                          node_color=cmap(norm((self.attribute))),
                                          alpha=1,
                                          node_size = [size * n for n in nsize])

        nodes.set_edgecolor('grey')
        nodes.set_linewidth(2)
        nx.draw_networkx_edges(self.graph,
                               pos = pos,
                               width= 3,
                               alpha = 0.5,
                               edge_color='dimgray')


       #position labels slighter offset to nodes
        if labels == True:
            pos_higher = {}
            off = 0.05  # offset on the y axis

            for k, v in pos.items():
                pos_higher[k] = (v[0], v[1])

            nx.draw_networkx_labels(self.graph,
                                    pos_higher,
                                    font_size=20,
                                    font_color="black",
                                    bbox = {"ec": "k", "fc": "white", "alpha": 0.6})


        #legend
        ticks = [min(self.attribute), max(self.attribute)]
        ax = plt.gca()
        ax.axis('off')
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        cax = ax.inset_axes([1.1, 0.05, 0.04, 0.2])
        cbar = plt.colorbar(sm,cax=cax, ticks=ticks)
        cbar.ax.set_title(col_legend_title, fontsize = 22, pad = 30)
        cbar.ax.tick_params(labelsize=22, pad =10, bottom = True)

        #define the size of the node labels and plot as non-existent dots
        #obtain the 10th, 50th, and 100th percentile
        legend_n = {int(round(np.percentile(list(self.samplecount_node.values()), 10),0)),
                    int(round(np.percentile(list(self.samplecount_node.values()), 50),-1)),
                    int(round(np.percentile(list(self.samplecount_node.values()), 100),-1))}

        for v in legend_n:
            plt.scatter([],[], s= (size * v), label='{}'.format(v))

        #get the legend handles and arrange to correct order
        handles,labels = ax.get_legend_handles_labels()
        handles, labels = zip(*[ (handles[i], labels[i]) for i in sorted(range(len(handles)), key=lambda k: list(map(int,labels))[k])] )
        ax.legend(handles,labels,loc=2, bbox_to_anchor=(1, 1),fancybox=True, title = "# of samples", labelspacing  = 1.6,frameon=False)

        #set colour of points
        leg = ax.get_legend()
        for i in [0,1,2]:
            try:
                leg.legendHandles[i].set_color('grey')
            except Exception:
                pass

        plt.show()

        self.node_colour_legend = cmap(norm((self.attribute)))
