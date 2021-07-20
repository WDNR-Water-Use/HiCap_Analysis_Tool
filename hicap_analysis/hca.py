from wells import Well, WellResponse
import numpy as np
import pandas as pd

class Hca():
    def __init__(self, T, S, time=[], Q=[],stream_locs={},  stream_apportionment={}, assessed_well_locs=[]) -> None:
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
            assessed_well_locs (list, optional): [description]. Defaults to [].
        """
        pass

    def existing_wells(self):
        """dictionary of existing Well objects for cumulative impacts
        """

        pass

    def proposed_wells(self):
        """dictionary of proposed Well object(s)
        """
        pass

    def existing_impacts(self):
        """possibly precalculated so a load or calculation option
        """
        
        pass

    def proposed_impacts(self):
        """
        calculate and store all the new impacts
        """
        pass

    def calculate_significance(self):
        """[summary]
        """
        pass


