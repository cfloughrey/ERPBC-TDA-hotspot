import hotmapper.mapper as hm
import hotmapper.utils as hmu
import hotmapper.hotspot as hmh
import hotmapper.random_lens as hmrl

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
        self.runs = runs
        self.parameters = {}
        self.parameter_lens = []
        self.parameter_samples = {}

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
                random_lens = hmrl.linear(self.X, n_features = parameters["lens_feature_no"])

                #build a grid of interval and overlap combinations
                io_list = list(product(parameters["interval_list"], parameters["overlap_list"]))

                #for each grid combination of the interval & overlap, search lens for hotspot
                for i in io_list:
                    #define intervals and overlap
                    i_param = i[0]
                    o_param = i[1]

                    # Run a clustering algorithm and build the graph
                    mapper = hm.Mapper(data = self.X,
                                            lens_function = random_lens["lens"],
                                            intervals = i_param,
                                            overlap = o_param,
                                            clustering_algorithm = parameters["clustering_algorithm"],
                                            text = False)

                    #build the graph with edges and nodes
                    mapper.build_graph()

                    #visualise graph
                    if visualise == True:
                        hmv.draw_graph(mapper_graph = mapper.graph,
                                        attribute_function = parameters["attribute_function"],
                                        samples_in_nodes = mapper.samples_in_nodes,
                                        size = 5,
                                        style = 2,
                                        labels = False)

                    #run hotspot detection
                    hotspot_search = hmh.HotSpot(mapper_graph = mapper.graph,
                                                 attribute_function = parameters["attribute_function"],
                                                 samples_in_nodes = mapper.samples_in_nodes)

                    hotspots = hotspot_search.search_graph(attribute_threshold = parameters["epsilon"],
                                                                min_sample_size = parameters["min_samples"],
                                                                attribute_extreme = parameters["extreme"])

                    #return list of samples in each hotspot found
                    sample_list = []
                    for n in hotspots:
                        sample_list.append(hmu.sample_index_in_nodes(mapper.samples_in_nodes, n))


                    #if hotspot present, save properties
                    if any(hotspots):
                        self.parameters[(i_param,o_param)] = hotspots # list of hotspots
                        self.parameter_lens = [random_lens["weights"],random_lens["feature_list"]]
                        self.parameter_samples[(i_param,o_param)] = sample_list
                        significance = True

                #if hotspots exist in the filter function search
                if self.parameters:
                    print("\nHotspots search successful")
                    return
                else:
                    count += 1
                    print(F"\nCompleted {count} searches")
