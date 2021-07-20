import numpy as np
from numpy.lib.arraysetops import isin
import pandas as pd
import scipy.special as sps

# define drawdown methods here
def _theis(T,S,time,dist,Q):
    """Calculate Theis drawdown at single location

    Args:
    
    T (float): transmissivity [ft**2/d]
    S (float): storage [unitless]
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float, optionall np.array or list): distance at which to calculate results in [ft]
    Q (float): pumping rate (+ is extraction) [ft**3/d]
    
    Returns:
        float (array): drawdown values at input parameter
                        times/distances 
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
    """Calcualte Glover 

    Args:
    T (float): transmissivity [ft**2/d]
    S (float): storage [unitless]
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float, optionall np.array or list): distance at which to calculate results in [ft]
    Q (float): pumping rate (+ is extraction) [ft**3/d]


    Returns:
        float (array): depletion values at at input parameter
                        times/distances
    """
    z = np.sqrt((dist**2*S)/(4 * T * time))
    return Q * sps.erfc(z)

def _sdf(T,S,dist,**kwargs):
    """[summary]

    Args:
        T (float): transmissivity [ft**2/d]
        S (float): storage [unitless]
        dist (float, optionall np.array or list): distance at which to calculate results in [ft]
        **kwargs: just included to all for extra values in call
    Returns:
        SDF: Stream depletion factor [d]
    """
    if isinstance(dist, list):
        dist = np.array(dist)
    return dist**2*S/T

def _walton(T,S,dist,time, Q):
    """Calcualte depletion using Watkins (1987) PT-8 BASIC program logic 

    Args:
    T (float): transmissivity [gpd/ft]
    S (float): storage [unitless]
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float): distance at which to calculate results in [ft]
    Q (float): pumping rate (+ is extraction) [ft**3/d]


    Returns:
        float (array): depletion values at at input parameter
                        times/distances [CFS]
    """
    if dist > 0:
        G = dist / np.sqrt((0.535 * time * T/S))
    else:
        G = 0
    I = 1 + .0705230784*G + .0422820123*(G**2) + 9.2705272e-03*(G**3)
    J = (I + 1.52014E-04*(G**4) + 2.76567E-04*(G**5)+4.30638E-05*(G**6)) ** 16
    return Q * (1/J) / 3600 / 24

ALL_DD_METHODS = {'theis': _theis}

ALL_DEPL_METHODS = {'glover': _glover,
                    'walton': _walton}
GPM2CFD = 192.5

class WellResponse():
    """[summary]
    """
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
    """[summary]
    """
    '''
    object to contain a proposed or existing well with all responses for that one well
    defined and to preprocess
    '''
    def __init__(self, well_loc, T, S, Q, time, stream_locs={},  stream_apportionment={}, assessed_well_locs=[]) -> None:
        """[summary]

        Args:
            well_loc ([type]): [description]
            T ([type]): [description]
            S ([type]): [description]
            Q ([type]): [description]
            time ([type]): [description]
            stream_locs (dict, optional): [description]. Defaults to {}.
            stream_apportionment (dict, optional): [description]. Defaults to {}.
            assessed_well_locs (list, optional): [description]. Defaults to [].
        """
        self.well_loc=well_loc
        self.stream_locs=stream_locs
        self.assessed_well_locs=assessed_well_locs # wells at which defining impacts (drawdown)
        self.T = T
        self.S = S
        self.time = time
        self.Q = Q
        self.stream_apportionment=stream_apportionment
    
    def preprocess(self):
        """        
        -convert locations to distances
        -instantiate a WellResponse object for each response (each stream and each assessed well)
        """

        pass

    def calc_responses(self):
        """populate the depletion and drawdown proprties (calc if needed)
        """
        
        
        self.depletion = {} # keys are response names, values are arrays of depletion
        self.drawdown = {} # keys are response names, values are arrays of drawdown
        
