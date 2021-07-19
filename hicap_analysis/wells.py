import numpy as np
from numpy.lib.arraysetops import isin
import pandas as pd
import scipy.special as sps

# define drawdown methods here
def _theis(T,S,time,dist,Q):
    """
    T (float): transmissivity [ft**2/d]
    S (float): storage [unitless]
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float, optionall np.array or list): distance at which to calculate results in [ft]
    Q (float): pumping rate (+ is extraction) [ft**3/d]
    """
    if isinstance(time, list):
        time = np.array(time)
    if isinstance(dist, list):
        dist = np.array(dist)
    # contruct the well function argument
    u = dist**2. * S / (4. * T * time)
    # calculate and return
    return (Q / (4. * np.pi * T)) * sps.exp1(u)
    
# define stream depletion methods here
def _glover(T,S,time,dist,Q):
    #maths
    return None

ALL_DD_METHODS = {'theis': _theis}

ALL_DEPL_METHODS = {'glover': _glover}
GPM2CFD = 192.5

class WellResponse():
    '''
    Single well to calculate single response
    '''
    def __init__(self, name, T, S, dist, time, Q, stream_apportionment, dd_method='Theis', depl_method= 'Glover') -> None:
        self.name = name # name of response (stream, lake, or assessed well) evaluated
        self.T = T
        self.S = S
        self.dist = dist
        self.time = time
        self.dd_method=dd_method
        self.depl_method=depl_method
        self.drawdown=None
        self.depletion=None
        self.Q = Q
        self.stream_apportionment = stream_apportionment
        
    def _calc_drawdown(self):
        dd_f = ALL_DD_METHODS[self.dd_method.lower()]

        self.drawdown = 99999
        pass
    def _calc_depletion(self):
        dd_f = ALL_DEPL_METHODS[self.depl_method.lower()]
        self.depletion = 999999
        pass
    
    @property
    def drawdown(self):
        if self.drawdown is None:
            return self._calc_drawdown
        else:
            return self.drawdown

    @property
    def depletion(self):
        if self.depletion is None:
            return self._calc_depletion
        else:
            return self.depletion


class Well():
    '''
    object to contain a proposed or existing well with all responses for that one well
    defined and to preprocess
    '''
    def __init__(self, well_loc, T, S, Q, time, stream_locs={},  stream_apportionment={}, assessed_well_locs=[]) -> None:
        '''
        stream_locs (dict): keys are names, values are list-like locations
        assessed_well_locs (dict): keys are names, values are list-like locations
        
        most of this gets passed on to calculations for wells
        '''
        self.well_loc=well_loc
        self.stream_locs=stream_locs
        self.assessed_well_locs=assessed_well_locs # wells at which defining impacts (drawdown)
        self.T = T
        self.S = S
        self.time = time
        self.Q = Q
        self.stream_apportionment=stream_apportionment
    
    def preprocess(self):
        '''
        -convert locations to distances
        -instantiate a WellResponse object for each response (each stream and each assessed well)
        '''
        pass

    def calc_responses(self):
        '''
        populate the depletion and drawdown proprties (calc if needed)
        '''
        self.depletion = {} # keys are response names, values are arrays of depletion
        self.drawdown = {} # keys are response names, values are arrays of drawdown
        
