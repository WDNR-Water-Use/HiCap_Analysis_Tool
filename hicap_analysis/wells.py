import numpy as np
import pandas as pd
from pandas.core import base
import scipy.special as sps
import sys

# define drawdown methods here
def _theis(T,S,time,dist,Q, **kwargs):
    """Calculate Theis drawdown at single location

    Args:
    
    T (float): transmissivity [ft**2/d]
    S (float): storage [unitless]
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float, optionally np.array or list): distance at which to calculate results in [ft]
    Q (float): pumping rate (+ is extraction) [ft**3/d]
    **kwargs: just included to all for extra values in call    
    Returns:
        float (array): drawdown values at input parameter
                        times/distances 
    """
    if isinstance(time, list):
        time = np.array(time)
    if isinstance(dist, list):
        dist = np.array(dist)
    # construct the well function argument
    u = dist**2. * S / (4. * T * time)
    # calculate and return
    return (Q / (4. * np.pi * T)) * sps.exp1(u)
    
# define stream depletion methods here
def _glover(T,S,time,dist,Q, **kwargs):
    """Calculate Glover 
    from Glover and Balmer (1954)
    Args:
    T (float): transmissivity [ft**2/d] (K*D in the original paper)
    S (float): storage [unitless] (V in the original paper)
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float, optionally np.array or list): distance at which to calculate results in [ft] (X1 in the paper)
    Q (float): pumping rate (+ is extraction) [ft**3/d]
    **kwargs: just included to all for extra values in call

    Returns:
        float (array): depletion values at at input parameter
                        times/distances
    """
    z = dist/np.sqrt(4 * (T/S) * time)
    return Q * sps.erfc(z) /3600 / 24 # from CFD back to CFS

def _sdf(T,S,dist,**kwargs):
    """[summary]

    Args:
        T (float): transmissivity [ft**2/d]
        S (float): storage [unitless]
        dist (float, optionally np.array or list): distance at which to calculate results in [ft]
        **kwargs: just included to all for extra values in call
    Returns:
        SDF: Stream depletion factor [d]
    """
    if isinstance(dist, list):
        dist = np.array(dist)
    return dist**2*S/T

def _walton(T,S,time,dist, Q,**kwargs):
    """Calculate depletion using Watkins (1987) PT-8 BASIC program logic 

    Args:
    T (float): transmissivity [gpd/ft]
    S (float): storage [unitless]
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float): distance at which to calculate results in [ft]
    Q (float): pumping rate (+ is extraction) [ft**3/d]
    **kwargs: just included to all for extra values in call

    Returns:
        float (array): depletion values at at input parameter
                        times/distances [CFS]
    """
    if isinstance(time, list):
        time = np.array(time)
    if isinstance(time, pd.Series):
        time=time.values
    if dist > 0:
        # avoid divide by zero for time==0
        # time = time.values
        G = np.zeros_like(time).astype(float)
        G[time!=0] = dist / np.sqrt((0.535 * time[time!=0] * T/S))
    else:
        G = 0
    I = 1 + .0705230784*G + .0422820123*(G**2) + 9.2705272e-03*(G**3)
    J = (I + 1.52014E-04*(G**4) + 2.76567E-04*(G**5)+4.30638E-05*(G**6)) ** 16
    ret_vals = Q * (1/J) / 3600 / 24 # from CFD back to CFS
    ret_vals[time==0] = 0.0
    return ret_vals

def _hunt99(T,S,time,dist,Q, **kwargs):
    '''Function for Hunt (1999) solution for streamflow depletion
    by a pumping well.  

    Hunt, B., 1999, Unsteady streamflow depletion from ground
    water pumping: Groundwater, v. 37, no. 1, pgs. 98-102, 
    https://doi.org/10.1111/j.1745-6584.1999.tb00962.x

    Parameters
    ----------
    T: float
        Transmissivity of aquifer (ft^2/day)
    S: float
        Storativity of aquifer (dimensionless)
    time (float, optionally np.array or list): time at which to calculate results [d]
    dist (float, optionally np.array or list): distance at which to calculate results in [ft] (l in the paper)
    use an array for either time or distance but not both
    Q (float): pumping rate (+ is extraction) [ft**3/d]
    streambed (float): streambed conductance [ft/d] (lambda in the paper)
    **kwargs: just included to all for extra values in call
    
    Returns
    -------
    Qs (float): streamflow depletion rate (CFS), optionally np.array or list 
                depending on input of time and dist
    '''
    if 'streambed' in kwargs.keys():
        streambed = kwargs['streambed']
    # turn lists into np.array so they get handled correctly
    if isinstance(time, list) and isinstance(dist, list):
        print('cannot have both time and distance as arrays')
        print('in the Hunt99 method.  Need to externally loop')
        print('over one of the arrays and pass the other')
        sys.exit()
    elif isinstance(time, list):
        time = np.array(time)
    elif isinstance(dist, list):
        dist = np.array(dist)

    a = np.sqrt(S * dist**2/(4. * T * time))
    b = (streambed**2 * time)/(4 * S * T)
    c = (streambed * dist)/(2. * T)
    # Qs/Q = erfc(a) - exp(b+c)*erfc(sqrt(b) + a)
    # in order to calculate exp(x)erfc(y) 
    # as values get big. Use numpy special erfcx(),
    # scaled complementary error function, which returns 
    # exp(y**2)erfc(y)
    # then compute exp(x)/exp(y**2) * erfcx(y)
    # which may be computed as exp(x - y**2) * erfcx(y)
    # This approach gives a product that goes to zero
    # as the exp() term gets big and erfc() goes to zero
    y = np.sqrt(b) + a
    t1 = sps.erfcx(y)
    t2 = np.exp(b+c-y**2)
    depl = sps.erfc(a) - (t1*t2)
    return (Q / (3600 * 24)) * depl

ALL_DD_METHODS = {'theis': _theis}

ALL_DEPL_METHODS = {'glover': _glover,
                    'walton': _walton,
                    'hunt99': _hunt99}

GPM2CFD = 60*24/7.48 # factor to convert from GPM to CFD

class WellResponse():
    """[summary]
    """
    def __init__(self, name, response_type, T, S, dist, Q, stream_apportionment=None, 
                    dd_method='Theis', depl_method= 'Glover', theis_time = -9999,
                    depl_pump_time = -99999, streambed_conductance=None) -> None:
        """Class to calculate a single response for a single pumping well.
        Args:
            name (str): pumping well name
            response_type (str): reserved for future implementation
            T (float): Aquifer Transmissivity
            S (float): Aquifer Storage
            dist (float): Distance between well and response
            Q (pandas Series): ...
            stream_apportionment ([type], optional): [description]. Defaults to None.
            dd_method (str, optional): [description]. Defaults to 'Theis'.
            depl_method (str, optional): [description]. Defaults to 'Glover'.
            theis_time (int, optional): [description]. Defaults to -9999.
            depl_pump_time (int, optional): [description]. Defaults to -99999.
            streambed_conductance (float, optional): Streambed conductance for the Hunt99 depletion
                        method. Defaults to None
        """
        self._drawdown=None
        self._depletion=None
        self.name = name # name of response (stream, or drawdown response (e.g. assessed well, lake, spring)) evaluated
        self.response_type = response_type # might use this later to sort out which response to return
        self.T = T
        self.T_gpd_ft = T*7.48
        self.S = S
        self.dist = dist
        self.dd_method=dd_method
        self.depl_method=depl_method
        self.theis_time = theis_time
        self.depl_pump_time = depl_pump_time
        self.Q = Q
        self.stream_apportionment = stream_apportionment
        self.streambed_conductance = streambed_conductance
        
    def _calc_drawdown(self):
        """calculate drawdown at requested distance and time
        """
        dd_f = ALL_DD_METHODS[self.dd_method.lower()]

        return dd_f(self.T, self.S, self.theis_time, self.dist, self.Q.iloc[0])
        
    def _calc_depletion(self):
        depl_f = ALL_DEPL_METHODS[self.depl_method.lower()]
        
        # find the differences in pumping
        dq = self.Q.copy()
        dq.iloc[1:] = np.diff(self.Q)
        
        # start with zero depletion
        depl = np.zeros(len(self.Q))
        
        # get the locations of changes
        deltaQ = dq.loc[dq!=0]
        
        # special case for starting with 0 pumping
        if self.Q.index[0] not in deltaQ.index:
            deltaQ.loc[self.Q.index[0]] = self.Q.iloc[0]
            deltaQ.sort_index(inplace=True)
        # initialize with pumping at the first time being positive
        idx = deltaQ.index[0]-1
        cQ = deltaQ.iloc[0]
        ct = list(range(idx,len(self.Q)))
        if self.depl_method.lower() == 'walton':
            # Walton method (only) needs these goofy units of gpd/dt for T
            T = self.T_gpd_ft
        else:
            T = self.T
        if self.depl_method.lower() == 'hunt99':
            extra_args = {'streambed':self.streambed_conductance}
        else:
            extra_args = {}
        depl[idx:] = depl_f(T, self.S, ct, self.dist, cQ*self.stream_apportionment, **extra_args)
        if len(deltaQ) > 1:
            deltaQ = deltaQ.iloc[1:]
            for idx,cQ in zip(deltaQ.index,deltaQ.values):
                idx-=2
                ct = list(range(len(self.Q)-idx))
                # note that by setting Q negative from the diff calculations, we always add
                # below for the image wells
                depl[idx:] += depl_f(T, self.S, ct, self.dist, cQ*self.stream_apportionment, **extra_args)
        return depl

    
    @property
    def drawdown(self):
        return self._calc_drawdown()

    @property
    def depletion(self):
        return self._calc_depletion()     

class Well():
    """Object to evaluate a pending (or existing, or a couple other possibilities) well with all relevant impacts.
        Preproceessing makes unit conversions and calculates distances as needed
    """

    def __init__(self, well_status='pending', T=-9999, S=-99999, Q=-99999, depletion_years=5, theis_dd_days=-9999, depl_pump_time=-9999,
         stream_dist=None, drawdown_dist=None,  stream_apportionment=None, depl_method='walton', streambed_conductance = None) -> None:
        """[summary]

        Args:
            T ([type]): [description]
            S ([type]): [description]
            Q ([type]): [description]
            depletion_years (int, optional): [description]. Defaults to 4.
            theis_dd_days (int, optional): [description]. Defaults to -9999.
            depl_pump_time (int, optional): [description]. Defaults to -9999.
            stream_dist ([type], optional): [description]. Defaults to None.
            drawdown_dist ([type], optional): [description]. Defaults to None.
            stream_apportionment ([type], optional): [description]. Defaults to None.
            depl_method ([str], optional): [description]. Defaults to walton
            streambed_conductance(dict: optional): dictionary of streambed conductances. Defaults to None
        """

        # placeholders for values returned with @property decorators
        self._depletion = None
        self._drawdown = None
        self._max_depletion = None
        self.depl_method = depl_method
        self.stream_dist = stream_dist
        self.drawdown_dist = drawdown_dist
        self.T = T
        self.S = S
        self.depletion_years = depletion_years
        self.theis_dd_days = theis_dd_days
        self.depl_pump_time = depl_pump_time
        self.Q = Q
        self.stream_apportionment=stream_apportionment
        self.streambed_conductance = streambed_conductance
        self.stream_responses = {} # dict of WellResponse objects for this well with streams
        self.drawdown_responses = {} # dict of WellResponse objects for this well with drawdown responses
        self.well_status = well_status # this is for the well object - later used for aggregation and must be
                # {'existing', 'active', 'pending', 'new_approved', 'inactive' }
        # make sure stream names consistent between dist and apportionment
        if stream_dist is not None and stream_apportionment is not None:
            assert len(set(self.stream_dist.keys())-set(self.stream_apportionment.keys())) == 0
        if stream_dist is not None and streambed_conductance is not None:
            assert len(set(self.stream_dist.keys())-set(self.streambed_conductance.keys())) == 0
        if self.stream_dist is not None:
            self.stream_response_names = list(self.stream_dist.keys())
        if self.drawdown_dist is not None:
            self.drawdown_response_names = list(self.drawdown_dist.keys())

        
        # now make all the WellResponse objects
        # first for streams
        if self.stream_dist is not None:
            for cs, (cname, cdist) in enumerate(self.stream_dist.items()):
                if self.streambed_conductance is not None:
                    streambed_conductance_current = self.streambed_conductance[cname]
                else:
                    streambed_conductance_current = None
                self.stream_responses[cs+1] = WellResponse(cname, 'stream', T=self.T, S=self.S, 
                                    dist=cdist, depl_pump_time =self.depl_pump_time, 
                                    Q=self.Q, stream_apportionment=self.stream_apportionment[cname], 
                                    depl_method=self.depl_method, 
                                    streambed_conductance=streambed_conductance_current)

        # next for drawdown responses
        if self.drawdown_dist is not None:
            for cw, (cname,cdist) in enumerate(self.drawdown_dist.items()):
                self.drawdown_responses[cw+1] = WellResponse(cname, 'well', T=self.T, S=self.S, 
                                    dist=cdist, theis_time=self.theis_dd_days, 
                                    Q=self.Q, dd_method='theis')      
        
    @property
    def drawdown(self):
        if self._drawdown is None:
            self._drawdown = {}
            for cw, cwob in self.drawdown_responses.items():
                self._drawdown[cwob.name] = cwob.drawdown
        return self._drawdown

    @property
    def depletion(self):
        self._depletion = {}
        for _ , cwob in self.stream_responses.items():
            self._depletion[cwob.name] = cwob.depletion
           
        return self._depletion

    @property
    def max_depletion(self):
        return {cwob.name:np.max(cwob.depletion) 
                for _,cwob in self.stream_responses.items()}
        

    

