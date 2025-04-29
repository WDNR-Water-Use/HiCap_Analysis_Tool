import numpy as np
import pandas as pd
from pandas.core import base
import scipy.special as sps
import scipy.integrate as integrate
from scipy.special import gammaln
import sys
import hicap_analysis

class WellResponse():
    """Class to facilitate depletion or drawdown calculations

    Objects from this class will have the required information needed to
    call any of the analysis methods in the package.  The intent is that
    users will not have to access analytical solutions directly, but 
    can set up a WellResponse object.  The generation of WellResonse
    objects is generally done through an AnalysisProject object.

    """
    def __init__(self, name, response_type, T, S, dist, Q, stream_apportionment=None, 
                    dd_method='Theis', depl_method= 'Glover', theis_time = -9999,
                    depl_pump_time = -99999, streambed_conductance=None) -> None:
        """Class to calculate a single response for a single pumping well.

        *refactor for kwargs?*

        Parameters
        ----------
        name: string 
            pumping well name
        response_type: string
            reserved for future implementation
        T: float
            Aquifer Transmissivity
        S: float
            Aquifer Storage
        dist: float
            Distance between well and response
        Q: pandas series
            Pumping rate changes and times
        stream_apportionment: string
                ([type], optional): [description]. Defaults to None.
        dd_method: string
            [description]. Defaults to 'Theis'.
        depl_method:string 
            (str, optional): [description]. Defaults to 'Glover'.
        theis_time: integer
            (int, optional): [description]. Defaults to -9999.
        depl_pump_time: integer
            (int, optional): [description]. Defaults to -99999.
        streambed_conductance: float
            Streambed conductance for the Hunt99 depletion method. Defaults to None

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
        """calculate drawdown at requested distance and time using solution given as attribute to the object
        """
        dd_f = hicap_analysis.ALL_DD_METHODS[self.dd_method.lower()]
        # start with zero drawdown
        dd = np.zeros(len(self.Q))
        deltaQ = hicap_analysis._calc_deltaQ(self.Q.copy())
        # initialize with pumping at the first time being positive
        idx = deltaQ.index[0]-1
        cQ = deltaQ.iloc[0]
        ct = list(range(idx,len(self.Q)))       
        extra_args= {}
        dd[idx:] = dd_f(self.T, self.S, ct, self.dist, cQ, **extra_args)
        if len(deltaQ) > 1:
            deltaQ = deltaQ.iloc[1:]
            for idx,cQ in zip(deltaQ.index,deltaQ.values):
                idx-=2
                ct = list(range(len(self.Q)-idx))
                # note that by setting Q negative from the diff calculations, we always add
                # below for the image wells
                dd[idx:] += dd_f(self.T, self.S, ct, self.dist, cQ, **extra_args)
        return dd
        
    def _calc_depletion(self):
        """calculate streamflow depletion at time using solution given as attribute to the object
        """
        depl_f = hicap_analysis.ALL_DEPL_METHODS[self.depl_method.lower()]
        # start with zero depletion
        depl = np.zeros(len(self.Q))
        
        deltaQ = hicap_analysis._calc_deltaQ(self.Q.copy())

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

        Parameters
        ----------
        T: float
            Aquifer Transmissivity
        S: float
            Aquifer Storage
        Q: pandas series
            Pumping rate changes and times
        depletion_years: int
           [description]. Defaults to 4.
        theis_dd_days: int
            [description]. Defaults to -9999.
        depl_pump_time: int 
            [description]. Defaults to -9999.
        stream_dist: type
            [description]. Defaults to None.
        drawdown_dist: type
            [description]. Defaults to None.
        stream_apportionment: type
            [description]. Defaults to None.
        depl_method: string
            description]. Defaults to walton
        streambed_conductance: dict
            dictionary of streambed conductances. Defaults to None
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
        

    

