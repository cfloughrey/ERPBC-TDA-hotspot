import hotmapper.mapper as hm
import hotmapper.utils as hmu
import hotmapper.hotspot as hmh
from itertools import product

class Parameter_Search():
    """This class searches across the Mapper parameters to identify the parameters that build a graph which contains a hotspot

    Parameters
    ----------

    X : pandas dataframe
        Dataset from which to obtain lens and build mapper graphs

    runs: into
        How many times to run the search for a lens
            """

    def __init__(self, X, runs = 50):
        self.X = X
        self.mapper = hm.Mapper(X, text = False)
        self.runs = runs
        self.parameter_hs = {}
        self.parameter_io = {}
        self.parameter_lens = []

    def build_graphs(self, parameters, visualise = False):
        """Search through the parameter options and build mapper graphs

        Parameters
        ----------

        parameters : dictionary
            Dictionary of parameter options.
            Must include lens / clustering algorithm / interval list / overlap list / epsilon / minimum sample size / attribute function / hotspot extremity
        """

        #Runs = lens space
        count = 0
        signficance = False

        while count < self.runs:
            if signficance == False:
                #generate a random lens from features
                random_lens = self.mapper.lens_function(selection = parameters["lens_option"])

                #build a grid of interval and overlap combinations
                io_list = list(product(parameters["interval_list"], parameters["overlap_list"]))

                #for each grid combination of the interval & overlap, search lens for hotspot
                for i in io_list:
                    #define intervals and overlap
                    i_param = i[0]
                    o_param = i[1]

                    # Build a cover on the lens function - specify the number of intervals and the percentage overlap
                    self.mapper.covering(intervals = i_param, overlap = o_param)

                    # Run a clustering algorithm and build the graph
                    self.mapper.cluster_data(algorithm = parameters["clustering_algorithm"])

                    #build the graph with edges and nodes
                    graph = self.mapper.build_graph(parameters["attribute_function"])

                    #visualise graph
                    if visualise == True:
                        self.mapper.visualise(size = 10, style = 2, labels = True)

                    #run hotspot detection
                    graph_components = hmh.Subgraphs(self.mapper, text = False)
                    hotspots = graph_components.run_hotspot_search(attribute_threshold = parameters["epsilon"],
                                                                    min_sample_size = parameters["min_samples"],
                                                                    extreme = parameters["extreme"])

                    #if hotspot present, save properties
                    for h in hotspots:
                        check = hotspots[h]['hotspot class']
                        if any(check):
                            self.parameter_hs[h] = hotspots # list of hotspots
                            self.parameter_io[h] = [i_param,o_param]
                            self.parameter_lens = [self.mapper.random_value, self.mapper.random_index]
                            significance = True

                #if hotspots exist in the filter function search
                if self.parameter_hs:
                    print("\nHotspots search successful")
                    return
                else:
                    count += 1
                    print(F"\nCompleted {count} searches")
