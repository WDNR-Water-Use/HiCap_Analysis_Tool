from wells import Well, WellResponse
import numpy as np
import pandas as pd
def loc_to_dist(loc0, loc1):
    TODO: maths
    dist = loc0-loc1
    return dist
class Project():
    def __init__(self) -> None:
        """made up of multiple Well objects
        instantiate Well objects for existing and proposed
        have all names, properties, locations, etc.

        Args:
            T ([type]): [description]
            S ([type]): [description]
            time (list, optional): [description]. Defaults to [].
            Q (list, optional): [description]. Defaults to [].
            stream_locs (dict, optional): [description]. Defaults to {}.
            stream_apportionment (dict, optional): [description]. Defaults to {}.
            drawdown_locs (list, optional): [description]. Defaults to [].
        """
        pass

    def populate_from_yaml(self, ymlfile):
        pass
        # a bunch of Well objects

    def populate(self, T = -9999, S = -9999, existing_wells = {}, ):
        pass
        # a bunch of Well objects
        
    def aggregate_responses(self):
        # identify which responses and which wells eg. all of a certain status, or all,
        # or a single well
        pass



    def calculate_significance(self):
        """[summary]
        """
        pass


