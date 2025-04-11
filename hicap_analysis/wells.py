import numpy as np
import pandas as pd
from pandas.core import base
import scipy.special as sps
import scipy.integrate as integrate
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
    isarray = False
    if isinstance(time, list):
        time = np.array(time)
    if isinstance(dist, list):
        dist = np.array(dist)
    if isinstance(dist, pd.Series):
        dist = dist.values

    if isinstance(dist, np.ndarray):
        isarray = True
    
    # construct the well function argument
    # is dist is zero, then function does not exist
    # trap for dist==0 and set to small value
    if isarray:
        dist = np.where(dist==0, 0.001, dist)
    else:
        if dist == 0.:
            dist = 0.001
    
    # compute u and then the Theis solution
    u = dist**2. * S / (4. * T * time)
    
    # calculate and return
    return (Q / (4. * np.pi * T)) * sps.exp1(u)

def _hunt99ddwn(T, S, time, l, Q, **kwargs):
    ''' Calculate drawdown in an aquifer with a partially penetrating stream and
        including streambed resistance (Hunt, 1999).  The solution becomes the Theis
        solution if streambed conductance is zero, and approaches an image-well solution
        from Theis or Glover and Balmer (1954) as streambed conductance gets very large.
        Note that the well is located at the location x,y = (dist, 0) and the stream
        is aligned with y-axis at x=0 
    
    Parameters
    ----------
    T: float
        Transmissivity of aquifer (ft^2/day)
    S: float
        Storativity of aquifer (dimensionless)
    time: (float, optionally np.array or list): time at which to calculate results [d]
    l: distance between well and stream in [ft] 
    Q (float): pumping rate (+ is extraction) [ft**3/d]
    streambed (float): streambed conductance [ft/d] (lambda in the paper)
    x, y: either a pair of single x, y values to compute results or 
            vectors from numpy meshgrid giving grid of x,y locations
    **kwargs: just included to all for extra values in call
    
    Returns
    -------
    drawdown, or meshgrid of drawdowns, or np.array with shape (ntimes, meshgridxx, meshgridyy)

    '''
    if 'streambed' in kwargs.keys():
        streambed = kwargs['streambed']
    if 'x' in kwargs.keys():
        x = kwargs['x']
    if 'y' in kwargs.keys():
        y = kwargs['y']
        
    # turn lists into np.array so they get handled correctly,
    # check if time or space is an array
    timescalar = True
    spacescalar = True
    if isinstance(time, list):
        time = np.array(time)
    
    if isinstance(time, np.ndarray):
        timescalar = False

    if isinstance(x, np.ndarray):
        spacescalar = False
    
    # compute a single x, y point at a given time
    if timescalar and spacescalar:
        [strmintegral, err] = integrate.quad(_ddwn2, 0.0, np.inf,
                                         args=(l, x, y, T, streambed, time, S))
        return (Q / (4. * np.pi * T)) * (_ddwn1(l, x, y, T, streambed, time, S) - strmintegral)
    
    # compute a vector of times for a given point
    if not timescalar and spacescalar:
        drawdowns = []
        for tm in time:
            [strmintegral, err] = integrate.quad(_ddwn2, 0.0, np.inf,
                                         args=(l, x, y, T, streambed, tm, S))
            drawdowns.append((Q / (4. * np.pi * T)) * (_ddwn1(l, x, y, T, streambed, tm, S) - strmintegral))
        return drawdowns
    
    # if meshgrid is passed, return an np.array with dimensions
    # ntimes, num_x, num_y
    if not spacescalar:
        numrow = np.shape(x)[0]
        numcol = np.shape(x)[1]
        if timescalar:
            time = np.array([time])
        drawdowns = np.zeros(shape=(len(time), numrow, numcol))
        for time_idx in range(0, len(time)):
            for i in range(0, numrow):
                for j in range(0, numcol):
                    [strmintegral, err] = integrate.quad(_ddwn2, 0.0, np.inf,
                                         args=(l, x[i,j], y[i,j], T, streambed, time[time_idx], S))
                    drawdowns[time_idx, i, j] = (Q / (4. * np.pi * T)) * (_ddwn1(l, x[i,j], y[i,j], T, streambed, time[time_idx], S) - strmintegral)
        return drawdowns
    
def _ddwn1(l, x, y, T, streambed, time, S):
    ''' calculates Theis drawdown function for a point (x,y) given
        a well at the location (l, 0) from a stream.  Used in 
        computing Hunt, 1999 estimate of drawdown.  Equation 30 from 
        the paper.  Variables described in _hunt99ddwn function.
    '''
    if isinstance(l, list):
        l = np.array(l)
    if isinstance(l, pd.Series):
        l = l.values
    
    isarray = False
    if isinstance(l, np.ndarray):
        isarray = True
    
    # construct the well function argument
    # if (l-x) is zero, then function does not exist
    # trap for (l-x)==0 and set to small value
    dist = l - x
    if isarray:
        dist = np.where(dist==0, 0.001, dist)
    else:
        if dist == 0.:
            dist = 0.001

    u1 = ((dist)**2 + y**2)/(4. * T * time/S)
    
    return sps.exp1(u1)


def _ddwn2(theta, l, x, y, T, streambed, time, S):
    ''' calculates function that gets integrated in the Hunt (1999) drawdown
        equation (Equation 29 and 30 in the paper), theta is the constant
        of integration and the rest of the variables described in the
        _hunt99ddwn function.
    '''
    if streambed == 0.:
        return 0.
    u2 = ((l + np.abs(x) + 2*T*theta/streambed)**2 + y**2)/(4. * T * time/S)
    return np.exp(-theta) * sps.exp1(u2)

    
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

def _hunt99(T,S,time,dist,Q,streambed, **kwargs):
    ''' Function for Hunt (1999) solution for streamflow depletion
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
    dist (float, optionally np.array or list): distance between well and stream in [ft] (l in the paper)
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

def _hunt2003(T,S,time,dist,Q,Bprime, Bdouble, K, sigma, width, streambed, **kwargs):
    ''' Function for Hunt (2003) solution for streamflow depletion
    by a pumping well in a semiconfined aquifer.

    Hunt, B., 2003, Unsteady streamflow depletion when pumping
    from semiconfined aquifer: Journal of Hydrologic Engineering,
    v.8, no. 1, pgs 12-19. https://doi.org/10.1061/(ASCE)1084-0699(2003)8:1(12)
    

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
    Bprime: float
        saturated thickness of semiconfining layer containing stream, [ft]
    Bdouble: float
        distance from bottom of stream to bottom of semiconfining layer, [ft] (aquitard thickness beneath the stream)
    K: float
        hydraulic conductivity of semiconfining layer [ft/day]
    sigma: float
        porosity of semiconfining layer
    width: float
        stream width (b in paper) ,[ft]
    streambed (float): streambed conductance [ft/d] (lambda in the paper),
                        only used if K is less than 1e-10
    **kwargs: just included to all for extra values in call

    Returns
    -------
    Qs (float): streamflow depletion rate (CFS), optionally np.array or list 
                depending on input of time and dist
    '''
    # turn lists into np.array so they get handled correctly
    if isinstance(time, list) and isinstance(dist, list):
        print('cannot have both time and distance as arrays')
        print('in the Hunt2003 method.  Need to externally loop')
        print('over one of the arrays and pass the other')
        sys.exit()
    elif isinstance(time, list):
        time = np.array(time)
    elif isinstance(dist, list):
        dist = np.array(dist)

    # make dimensionless group used in equations
    dtime = (T * time) / (S * np.power(dist, 2))

    # if K is really small, set streambed conductance to a value
    # so solution collapses to Hunt 1999 (confined aquifer solution)
    if K<1.e-10:
        lam = streambed
    else:
        lam = K * width/Bdouble
    dlam = lam * dist/T
    epsilon = S/sigma
    dK = (K/Bprime) * np.power(dist, 2)/T

    # numerical integration of F() and G() functions to
    # get correction to Hunt(1999) estimate of streamflow depletion
    # because of storage in the semiconfining aquifer
    correction = []
    for dt in dtime:
        [y, err] = integrate.quad(_integrand, 
                                                0., 
                                                1., 
                                                args=(dlam, dt, epsilon, dK),
                                                limit=500)
        correction.append(dlam * y)
    
    # terms for depletion, similar to Hunt (1999) but repeated
    # here so it matches the 2003 paper.
    a = (1. / (2. * np.sqrt(dtime)))
    b = (dlam/2. + (dtime * np.power(dlam,2)/4.))
    c = a + (dlam * np.sqrt(dtime)/2.)

    # use erfxc() function from scipy (see _hunt99 above)
    # for erf(b)*erfc(c) term
    t1 = sps.erfcx(c)
    t2 = np.exp(b-c**2)
    depl = sps.erfc(a) - (t1*t2) 

    ## corrected depletion for storage of upper semiconfining unit
    return (Q / (3600 * 24)) * (depl - correction)           


def _F(alpha, dlam, dtime):
    ''' F function from paper in equation (46) as given
        by equation (47)

        Parameters
        ----------
        alpha: float
            integration variable
        dlam: float
            dimensionless streambed/semiconfining unit conductance
            (width * K/B'') * distance/Transmissivity
        dt: float
            dimensionless time
            (time * transmissivity)/(storativity * distance**2)
    '''
    # Hunt uses an expansion if dimensionless time>3 
    z = alpha*dlam*np.sqrt(dtime)/2. + 1./(2.*alpha*np.sqrt(dtime))
    if np.abs(z) < 3.0:
        a = dlam/2. + (dtime * np.power(alpha,2) * np.power(dlam,2)/4.)
        t1 = sps.erfcx(z)
        t2 = np.exp(a-z**2)
        b = -1./(4 * dtime * alpha**2)
        # equation 47 in paper
        F = np.exp(b) * np.sqrt(dtime/np.pi) - (alpha * dtime * dlam)/2. * (t1*t2)
    else:
        t1 = np.exp(-(1./(4.*dtime*alpha**2)))/(2.*alpha*z*np.sqrt(np.pi))
        t2 = 2./(dlam*(1.+(1./(dlam*dtime*alpha**2))**2))
        sumterm = 1 - (3./(2 * z**2)) + (15./(4. * z**4)) - (105./(8 * z**6))
        F = t1*(1.0 + t2*sumterm)  # equation 53 in paper

    if np.isnan(F):
        print(f'alpha {alpha}, dtime {dtime}, dlam {dlam}')
        sys.exit()
    return F

def _G(alpha, epsilon, dK, dtime):
    ''' G function from paper in equation (46) as given
        by equation (53). Uses scipy special for 
        incomplete Gamma Function (P(a,b)), binomial coefficient, 
        and modified Bessel function of zero order (I0).

        Parameters
        ----------
        alpha: float
            integration variable
        epsilon: float
            dimensionless storage
            storativity/porosity of semiconfining bed
        dK: float
            dimensionless conductivity of semiconfining unit
            (K * Bprime) * dist**2/Transmissivity
    '''
    # if dimensionless K is zero (check really small), return 0
    # this avoids divide by zero error in terms that have divide by (a+b)
    if dK < 1.0e-10:
        return 0.
    
    a = epsilon * dK * dtime * (1. - alpha**2)
    b = dK * dtime * alpha**2

    if (a + b) < 80.:
        term1 = np.exp(-(a+b))*sps.i0(2.*np.sqrt(a*b))
    else:
        term1 = 0.
    abterm = np.sqrt(a*b)/(a+b)

    sum = 0 
    for n in range(0, 101):
        if n <=8:
            addterm = sps.binom(2*n, n)*sps.gammainc(2*n+1,a+b)*abterm**(2*n)
        else:
            bi_term = np.log(sps.binom(2*n, n))
            inc_gamma = np.log(sps.gammainc(2*n+1, a+b))
            logab = (2*n)*np.log(abterm)
            addterm = np.exp(bi_term + inc_gamma + logab)
        sum = sum + addterm
        if addterm < 1.0e-08:
            break
    
    eqn52= 0.5 * (1. - term1 + ((b-a)/(a+b)) * sum)
    if eqn52 < 0:
        eqn52 = 0.
    if eqn52 > 1.:
        eqn52 = 1.

    if np.isnan(eqn52):
        print('equation 52 is nan')
        sys.exit()
    return eqn52


def _integrand(alpha, dlam, dtime, epsilon, dK):
    ''' Product of F() and G() terms for numerical 
        integration in equation 48
    
    '''
    return _F(alpha, dlam, dtime) * _G(alpha, epsilon, dK, dtime)

def _calc_deltaQ(Q):   
    """parse the Q time series to find changes and their associated times

    Args:
        Q (pandas Series): time series of pumping
    returns:
        deltaQ (pandas Series): times and changes in Q over time
    """
        # find the differences in pumping
    dq = Q.copy()
    dq.iloc[1:] = np.diff(Q)
        # get the locations of changes
    deltaQ = dq.loc[dq!=0]
    # special case for starting with 0 pumping
    if Q.index[0] not in deltaQ.index:
        deltaQ.loc[Q.index[0]] = Q.iloc[0]
        deltaQ.sort_index(inplace=True)
    return deltaQ


ALL_DD_METHODS = {'theis': _theis,
                    'hunt99ddwn': _hunt99ddwn}

ALL_DEPL_METHODS = {'glover': _glover,
                    'walton': _walton,
                    'hunt99': _hunt99,
                    'hunt03': _hunt2003}

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
        # start with zero drawdown
        dd = np.zeros(len(self.Q))
        deltaQ = _calc_deltaQ(self.Q.copy())
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
        depl_f = ALL_DEPL_METHODS[self.depl_method.lower()]
        # start with zero depletion
        depl = np.zeros(len(self.Q))
        
        deltaQ = _calc_deltaQ(self.Q.copy())

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
        

    

