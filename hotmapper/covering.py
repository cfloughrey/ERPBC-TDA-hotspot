import numpy as np
import pandas as pd

# This is included as a seperate class to the Mapper class so it can contain more relevant information
#than a simple dictionary attribute of the mapper class
class Cover():
    """This class builds a cover on the lens and identifies the corresponding samples in the dataset of each
    overlapping interval.
    """

    def __init__(self, data, lens, intervals, overlap):

        self.data = data
        self.lens = lens
        self.intervals = intervals
        self.overlap = overlap

        self.intervals_range = {}
        self.samples_in_interval = {}
        self.data_in_interval = {}



    def build_cover(self):
        """Divide the lens into the intervals and identify the samples in each interval
            """

        #find the range of the lens function
        lens_min = np.amin(self.lens)
        lens_max = np.amax(self.lens)

        #calculate the size of each interval
        interval_length = (lens_max - lens_min) / (((self.intervals-1)  *  (1 - self.overlap)) + 1)

        #the number of samples in the dataset
        no_samples = np.array(range(self.data.shape[0]))

        #these are the fundamental properties of the intervals
        interval_samples = [False for i in range(self.intervals)]
        interval_minmax = {}
        interval_samples = {}
        interval_data = {}


        for i in range (0,self.intervals):
            #the starting point of each interval is the start of the filter function
            ai = lens_min + (i * interval_length) * (1 - self.overlap)
            bi = ai + interval_length

            #each dict item contains the min and max of each interval
            interval_minmax[i] = [ai,bi]

            #for each point, assign it to an interval if the function value
            #of this point lies between interval start and end points
            points_in_set = ((ai <= self.lens) & (self.lens <= bi))
            samples = no_samples[points_in_set]
            interval_samples[i] = np.unique(samples)

            #return each value in lens_values if it is found in that set
            points = self.data[points_in_set]
            interval_data[i] = points



        #assign attribute of cover to class
        self.intervals_range = interval_minmax
        self.samples_in_interval = interval_samples
        self.data_in_interval = interval_data
