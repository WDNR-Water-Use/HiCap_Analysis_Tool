from wells import Well, WellResponse
import numpy as np
import pandas as pd

class Hca():
    def __init__(self, T, S, time=[], Q=[],stream_locs={},  stream_apportionment={}, assessed_well_locs=[]) -> None:
        '''
        made up of multiple Well objects
        instantiate Well objects for existing and proposed
        have all names, properties, locations, etc.
        '''
        pass

    def existing_wells(self):
        '''
        dictionary of existing Well objects for cumulative impacts
        '''
        pass

    def proposed_wells(self):
        '''
        dictionary of proposed Well object(s)
        pass
        '''
        pass

    def existing_impacts(self):
        '''
        possibly precalculated so a load or calculation option
        '''
        pass

    def proposed_impacts(self):
        '''
        calculate and store all the new impacts
        '''
        pass

    def calculate_significance(self):
        pass


