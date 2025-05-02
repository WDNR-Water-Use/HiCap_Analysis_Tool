import numpy as np
import pandas as pd
import scipy.special as sps
import scipy.integrate as integrate
from scipy.special import gammaln
import sys

""" FIle of drawdown and stream depletion analytical solutions
    as part of the pycap suite.

"""
# define drawdown methods here
def theis(T,S,time,dist,Q, **kwargs):
    """Calculate Theis drawdown 
    
        Calculates the Theis drawdown solution at specified times
        or distances from a pumping well.  

        https://pubs.usgs.gov/publication/70198446

   Parameters
   ----------
    T: float
        transmissivity [ft**2/d]
    S: float
        storage [unitless]
    time: float, optionally np.array or list
        time at which to calculate results [d]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft]
    Q: float
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs: included to all drawdown methods for extra values reuired in some calls    

    Returns
    -------
    drawdown: float or array of floats 
        drawdown values at input parameter times/distances [ft]

    Other Parameters
    ----------------
    **kwargs: dict
        no keyword arguments are needed for the Theis (1935) drawdown solution
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

def hunt99ddwn(T, S, time, l, Q, **kwargs):
    '''Calculate drawdown in an aquifer with a partially penetrating stream including streambed resistance (Hunt, 1999).  
        
        The solution becomes the Theis
        solution if streambed conductance is zero, and approaches an image-well solution
        from Theis or Glover and Balmer (1954) as streambed conductance gets very large.
        Note that the well is located at the location x,y = (l, 0) and the stream
        is aligned with y-axis at x=0 

        Hunt, B., 1999, Unsteady streamflow depletion from ground
        water pumping: Groundwater, v. 37, no. 1, pgs. 98-102, 
        https://doi.org/10.1111/j.1745-6584.1999.tb00962.x

    
    Parameters
    ----------
    T: float
        Transmissivity of aquifer (ft^2/day)
    S: float
        Storativity of aquifer (dimensionless)
    time: float, optionally np.array or list 
        time at which to calculate results [d]
    l: float
        distance between well and stream in [ft] 
    Q : float
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs:  See Other Parameters below
    
    Returns
    -------
    drawdown: float
        single value, meshgrid of drawdowns, or np.array with shape (ntimes, meshgridxx, meshgridyy)
        depending on input form of x, y, and ntimes

    Other Parameters
    ----------------
    **kwargs: dict
        keyword arguments needed for the Hunt (1999) drawdown solution are:
    streambed: float
        streambed conductance [ft/d] (lambda in the paper)
    x: float, optionally vector from numpy meshgrid giving grid of x,y locations
        distance from stream, [ft]
    y: float, optionally vector from numpy meshgrid giving grid of x,y locations
        distance parallel to stream, well is located at y=0
    
    call signature is:
        _hunt99ddwn(T, S, time, l, Q, streambed=streambed_value, x=X, y=Y)

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
    '''Internal method to calculate Theis drawdown function for a point (x,y) 
        
        Used in computing Hunt, 1999 estimate of drawdown.  Equation 30 from 
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
    '''Internal method to calculate function that gets integrated in the Hunt (1999) solution

        Equations 29 and 30 in the paper, theta is the constant
        of integration and the rest of the variables described in the
        _hunt99ddwn function.
    '''
    if streambed == 0.:
        return 0.
    u2 = ((l + np.abs(x) + 2*T*theta/streambed)**2 + y**2)/(4. * T * time/S)
    return np.exp(-theta) * sps.exp1(u2)

def WardLoughDrawdown(T1, T2, S1, S2, width, Q, dist, streambed_thick, streambed_K, aquitard_thick, aquitard_K, t, x, y, NSteh1=2, NSteh2=2):
    '''Compute drawdown using Ward and Lough (2011) solution

    *needs to be refactored to use kwargs*

        Ward and Lough (2011) presented a solution for streamflow depletion
        by a pumping well in a layered aquifer system.  The stream 
        is in the upper aquifer, and the pumping well is in a lower
        aquifer that is separated from the upper aquifer by a 
        semi-confining aquitard layer. 

        Ward, N.D.,and Lough, H., 2011, Stream depletion from pumping a 
        semiconfined aquifer in a two-layer leaky aquifer system (techical note):
        Journal of Hydrologic Engineering ASCE, v. 16, no. 11, pgs. 955-959,
        https://doi.org/10.1061/(ASCE)HE.1943-5584.0000382.

    Parameters
    ----------
    T: float 
        transmissivity [ft**2/d] (K*D in the original paper)
    S: float
        storage [unitless] (V in the original paper)
    time: float, optionally np.array or list
        time at which to calculate results [d]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft] (X1 in the paper)
    Q: float 
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs: see Other Parameters below
    
    Returns
    -------
    ddwn float
        drawdown at specified location or locations [ft]

    Other Parameters
    ----------------
    **kwargs: dict
        keyword arguments needed for the Hunt (2003) depletion solution are:
    T2: float
        Transmissivity of
    S2: float
        Storativity of
    streambed_thick: float
        thickness of streambed
    streambed_K: float
        hydraulic conductivity of streambed, ft/day
    aquitard_thick: float
        thickness of intervening leaky aquitard, ft
    aquitard_K: float
        hydraulic conductivity of intervening leaky aquifer, ft/day
    x: float
        xxxxxx
    y: float
        yyyyyy
    NSteh1: int
        IIIIII
    NStehl2: int
        IIIII
    width: float
        stream width (b in paper) [ft]
    
    '''
    # first nondimensionalize all the parameters
    x,y,t,T1,S1,K,lambd = _WardLoughNonDimensionalize(T1, T2, S1, S2, 
                                                    width, Q, dist, 
                                                    streambed_thick, streambed_K, 
                                                    aquitard_thick, aquitard_K, t, x, y)

    # Initialize output arrays
    s1 = np.zeros_like(t)
    s2 = np.zeros_like(t)

    # Inverse Fourier transform
    for ii in range(len(t)):
        try:
            s1[ii] = _StehfestCoeff(1, NSteh1) * _if1(T1, S1, K, lambd, x, y, np.log(2) / t[ii])
            for jj in range(2, NSteh1 + 1):
                s1[ii] += _StehfestCoeff(jj, NSteh1) * _if1(T1, S1, K, lambd, x, y, jj * np.log(2) / t[ii])
            s1[ii] *= np.log(2) / t[ii]
        except OverflowError as e:
            print(f"Overflow error in s1 calculation at index {ii}: {e}")
            s1[ii] = np.nan  # Assign NaN if there's an overflow

        try:
            s2[ii] = _StehfestCoeff(1, NSteh2) * _if2(T1, S1, K, lambd, x, y, np.log(2) / t[ii])
            for jj in range(2, NSteh2 + 1):
                s2[ii] += _StehfestCoeff(jj, NSteh2) * _if2(T1, S1, K, lambd, x, y, jj * np.log(2) / t[ii])
            s2[ii] *= np.log(2) / t[ii]
        except OverflowError as e:
            print(f"Overflow error in s2 calculation at index {ii}: {e}")
            s2[ii] = np.nan  # Assign NaN if there's an overflow

    return s1*Q/T2, s2*Q/T2 # re-dimensionalize


# define stream depletion methods here
def glover(T,S,time,dist,Q, **kwargs):
    """Calculate Glover and Balmer (1954) solution for stream depletion

        Depletion solution for a well near a river where the river fully
        penetrates the aquifer and there is no streambed resistance.
        
        Glover, R.E. and Balmer, G.G., 1954, River depletion from pumping
        a well near a river, Eos Transactions of the American Geophysical Union,
        v. 35, no. 3, pg. 468-470, https://doi.org/10.1029/TR035i003p00468.

    Parameters
    ----------
    T: float 
        transmissivity [ft**2/d] (K*D in the original paper)
    S: float
        storage [unitless] (V in the original paper)
    time: float, optionally np.array or list
        time at which to calculate results [d]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft] (X1 in the paper)
    Q: float 
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs: included to all depletion methods for extra values reuired in some calls 

    Returns
    -------
    drawdown: float 
        depletion values at at input parameter times/distances

    Other Parameters
    ----------------
    **kwargs: dict
        no keyword arguments are needed for the Glover and Balmer (1954) depletion solution    
    """
    z = dist/np.sqrt(4 * (T/S) * time)
    return Q * sps.erfc(z) /3600 / 24 # from CFD back to CFS

def sdf(T,S,dist,**kwargs):
    """internal function for Stream Depletion Factor

        Stream Depletion Factor was defined by Jenkins (1968) and described 
        in Jenkins as the time when the volume of stream depletion is 
        28 percent of the net volume pumped from the well.  
        SDF = dist**2 * S/T.

        Jenkins, C.T., Computation of rate and volume of stream depletion
        by wells: U.S. Geological Survey Techniques of Water-Resources
        Investigations, Chapter D1, Book 4, https://pubs.usgs.gov/twri/twri4d1/.

    Parameters
    ----------
    T: float
        transmissivity [ft**2/d]
    S: float
        storage [unitless]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft]
    **kwargs: just included to all for extra values in call

    Returns
    -------
    SDF: float
        Stream depletion factor [d]
    """
    if isinstance(dist, list):
        dist = np.array(dist)
    return dist**2*S/T

def walton(T,S,time,dist, Q,**kwargs):
    """Calculate depletion using Walton (1987) PT-8 BASIC program logic 

        Provides the Glover and Balmer (Jenkins) solution.

        Walton, W.C., Groundwater Pumping Tests:  Lewis Publishers, Chelsea, 
        Michigan, 201 p.

    Parameters
    ----------
    T: float 
        transmissivity [ft**2/d] (K*D in the original paper)
    S: float
        storage [unitless] (V in the original paper)
    time: float, optionally np.array or list
        time at which to calculate results [d]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft] (X1 in the paper)
    Q: float 
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs: included to all depletion methods for extra values required in some calls 

    Returns
    -------
    drawdown: float 
        depletion values at at input parameter times/distances

    Other Parameters
    ----------------
    **kwargs: dict
        no keyword arguments are needed for Walton (1987) depletion solution 
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

def hunt99(T,S,time,dist,Q,streambed, **kwargs):
    '''Function for Hunt (1999) solution for streamflow depletion by a pumping well. 

        Computes streamflow depletion by a pumping well for a partially penetrating
        stream with streambed resistance. 

        Hunt, B., 1999, Unsteady streamflow depletion from ground
        water pumping: Groundwater, v. 37, no. 1, pgs. 98-102, 
        https://doi.org/10.1111/j.1745-6584.1999.tb00962.x

    Parameters
    ----------
    T: float 
        transmissivity [ft**2/d] (K*D in the original paper)
    S: float
        storage [unitless] (V in the original paper)
    time: float, optionally np.array or list
        time at which to calculate results [d]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft] (X1 in the paper)
    Q: float 
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs: see Other Parameters below
    
    Returns
    -------
    Qs: float
        streamflow depletion rate, optionally np.array or list depending on input of time and dist [CFS]

    Other Parameters
    ----------------
    **kwargs: dict
        keyword arguments needed for the Hunt (1999) depletion solution are:
    streambed: float
        streambed conductance [ft/d] (lambda in the paper)
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

def hunt2003(T,S,time,dist,Q,Bprime, Bdouble, K, sigma, width, streambed, **kwargs):
    '''Function for Hunt (2003) solution for streamflow depletion by a pumping well.

        Computes streamflow depletion by a pumping well in a semiconfined aquifer
        for a partially penetrating stream.  The stream is in an upper semi-confining
        aquifer and pumping is in a lower aquifer.

        Hunt, B., 2003, Unsteady streamflow depletion when pumping
        from semiconfined aquifer: Journal of Hydrologic Engineering,
        v.8, no. 1, pgs 12-19. https://doi.org/10.1061/(ASCE)1084-0699(2003)8:1(12)
    

    Parameters
    ----------
    T: float 
        transmissivity [ft**2/d] (K*D in the original paper)
    S: float
        storage [unitless] (V in the original paper)
    time: float, optionally np.array or list
        time at which to calculate results [d]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft] (X1 in the paper)
    Q: float 
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs: see Other Parameters below
    
    Returns
    -------
    Qs: float
        streamflow depletion rate, optionally np.array or list depending on input of time and dist [CFS]

    Other Parameters
    ----------------
    **kwargs: dict
        keyword arguments needed for the Hunt (2003) depletion solution are:
    Bprime: float
        saturated thickness of semiconfining layer containing stream, [ft]
    Bdouble: float
        distance from bottom of stream to bottom of semiconfining layer, [ft] (aquitard thickness beneath the stream)
    K: float
        hydraulic conductivity of semiconfining layer [ft/day]
    sigma: float
        porosity of semiconfining layer
    width: float
        stream width (b in paper) [ft]
    streambed: float
        streambed conductance [ft/d] (lambda in the paper), only used if K is less than 1e-10
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
    '''F function from paper in equation (46) as given by equation (47) in Hunt (2003)

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
    '''G function from paper in equation (46) in Hunt (2003)

        This function is in equation (46) and expanded in
        equation (53). Function uses scipy special for 
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
    '''internal function returning product of F() and G() terms for numerical integration 
    
    '''
    return _F(alpha, dlam, dtime) * _G(alpha, epsilon, dK, dtime)

def _calc_deltaQ(Q):   
    """internal function to parse the Q time series to find changes and their associated times

    Parameters
    ----------
    Q: pandas Series
        time series of pumping

    Returns
    -------
    deltaQ: pandas Series)
        times and changes in Q over time
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


def _WardLoughNonDimensionalize(T1, T2, S1, S2, width, Q, dist, streambed_thick, streambed_K, aquitard_thick, aquitard_K, t, x, y):
    '''Internal function to make non-dimensional groups for Ward and Lough solution
    
    '''
    x /= dist
    y /= dist
    t = t*T2/(S2 * (dist**2))
    T1 /= T2
    S1 /= S2
    K = ((aquitard_K/aquitard_thick)*(dist**2))/T2
    lambd = ((streambed_K * width) / streambed_thick) *dist / T2
    return x,y,t,T1,S1,K,lambd    
    

def WardLoughDepletion(T1, T2, S1, S2, width, Q, dist, streambed_thick, streambed_K, aquitard_thick, aquitard_K, t, x=0, y=0,  NSteh1=2, NStehl2=2):
    '''Compute streamflow depletion using Ward and Lough (2011) solution

    *needs to be refactored to use kwargs*

        Ward and Lough (2011) presented a solution for streamflow depletion
        by a pumping well in a layered aquifer system.  The stream 
        is in the upper aquifer, and the pumping well is in a lower
        aquifer that is separated from the upper aquifer by a 
        semi-confining aquitard layer. 

        Ward, N.D.,and Lough, H., 2011, Stream depletion from pumping a 
        semiconfined aquifer in a two-layer leaky aquifer system (techical note):
        Journal of Hydrologic Engineering ASCE, v. 16, no. 11, pgs. 955-959,
        https://doi.org/10.1061/(ASCE)HE.1943-5584.0000382.

    Parameters
    ----------
    T: float 
        transmissivity [ft**2/d] (K*D in the original paper)
    S: float
        storage [unitless] (V in the original paper)
    time: float, optionally np.array or list
        time at which to calculate results [d]
    dist: float, optionally np.array or list
        distance at which to calculate results in [ft] (X1 in the paper)
    Q: float 
        pumping rate (+ is extraction) [ft**3/d]
    **kwargs: see Other Parameters below
    
    Returns
    -------
    Qs: float
        streamflow depletion rate, optionally np.array or list depending on input of time and dist [CFS]

    Other Parameters
    ----------------
    **kwargs: dict
        keyword arguments needed for the Hunt (2003) depletion solution are:
    T2: float
        Transmissivity of
    S2: float
        Storativity of
    streambed_thick: float
        thickness of streambed
    streambed_K: float
        hydraulic conductivity of streambed, ft/day
    aquitard_thick: float
        thickness of intervening leaky aquitard, ft
    aquitard_K: float
        hydraulic conductivity of intervening leaky aquifer, ft/day
    x: float
        xxxxxx
    y: float
        yyyyyy
    NSteh1: int
        IIIIII
    NStehl2: int
        IIIII
    width: float
        stream width (b in paper) [ft]
    
    '''
    # first nondimensionalize all the parameters
    x,y,t,T1,S1,K,lambd = _WardLoughNonDimensionalize(T1, T2, S1, S2, 
                                                    width, Q, dist, 
                                                    streambed_thick, streambed_K, 
                                                    aquitard_thick, aquitard_K, t, x, y)

    # Initialize output arrays
    s1 = np.zeros_like(t)
    s2 = np.zeros_like(t)
    # Inverse Fourier transform
    DeltaQ = _StehfestCoeff(1, NSteh1) * _if1_dQ(T1, S1, K, lambd, np.log(2) / t, x, y)
    for jj in range(2, NSteh1 + 1):
        DeltaQ += _StehfestCoeff(jj, NSteh1) * _if1_dQ(T1, S1, K, lambd, jj * np.log(2) / t, x, y)
    DeltaQ = 2 * np.pi * lambd * DeltaQ * np.log(2) / t

    return DeltaQ*Q

def _if1_dQ(T1, S1, K, lambda_, p, x, y):
    '''Internal function for Ward and Lough (2011) solution'''
    return _kernel1(T1,S1,K,lambda_,0,0,p)+_kernel2(T1,S1,K,lambda_,0,0,p)

def _if1(T1, S1, K, lambd, x, y, p):
    '''Internal function for Ward and Lough (2011) solution'''
    G = lambda phi: 2 * (_kernel1(T1, S1, K, lambd, x, np.tan(phi), p) +
                          _kernel2(T1, S1, K, lambd, x, np.tan(phi), p)) * \
                          np.cos(np.tan(phi) * y) / np.cos(phi) ** 2

    s1InvFour, _ = integrate.quad(G, 0, np.pi / 2, epsrel=1e-1, epsabs=1e-1, limit=10000)
    return s1InvFour

def _if2(T1, S1, K, lambd, x, y, p):
    '''Internal function for Ward and Lough (2011) solution'''
    H = lambda phi: 2 * (_coeff_s1_1(T1, S1, K, lambd, np.tan(phi), p) * 
                        _kernel1(T1, S1, K, lambd, x, np.tan(phi), p) +
                        _coeff_s1_2(T1, S1, K, lambd, np.tan(phi), p) *
                        _kernel2(T1, S1, K, lambd, x, np.tan(phi), p)) * \
                        np.cos(np.tan(phi) * y) / np.cos(phi) ** 2

    s2InvFour, errbnd = integrate.quad(H, 0, np.pi / 2, epsrel=1e-1, epsabs=1e-1, limit=10000)
    return s2InvFour

def _coeff_s1_1(T1, S1, K, lambd, theta, p):
    '''Internal function for Ward and Lough (2011) solution'''
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = _coeffs(T1, S1, K, lambd, theta, p)
    B1 = (mu1 * T1 - b11) / b12
    return B1

def _coeff_s1_2(T1, S1, K, lambd, theta, p):
    '''Internal function for Ward and Lough (2011) solution'''
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = _coeffs(T1, S1, K, lambd, theta, p)
    B2 = (mu2 * T1 - b11) / b12
    return B2

def _kernel1(T1, S1, K, lambd, x, theta_or_y, p):
    '''Internal function for Ward and Lough (2011) solution'''
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = _coeffs(T1, S1, K, lambd, theta_or_y, p)

    if x < 0:
        F1 = A1 * np.exp(x * np.sqrt(mu1))
    elif 0 <= x <= 1:
        F1 = A1 * np.exp(-x * np.sqrt(mu1)) + beta1 / (2 * np.sqrt(mu1) * l1) * \
             (np.exp((x - 1) * np.sqrt(mu1)) - np.exp(-(x + 1) * np.sqrt(mu1)))
    else:
        F1 = A1 * np.exp(-x * np.sqrt(mu1)) + beta1 / (2 * np.sqrt(mu1) * l1) * \
             (np.exp((1 - x) * np.sqrt(mu1)) - np.exp(-(x + 1) * np.sqrt(mu1)))
    return F1

def _kernel2(T1, S1, K, lambd, x, theta_or_y, p):
    '''Internal function for Ward and Lough (2011) solution'''
    b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2 = _coeffs(T1, S1, K, lambd, theta_or_y, p)

    if x < 0:
        F2 = A2 * np.exp(x * np.sqrt(mu2))
    elif 0 <= x <= 1:
        F2 = A2 * np.exp(-x * np.sqrt(mu2)) + beta2 / (2 * np.sqrt(mu2) * l2) * \
             (np.exp((x - 1) * np.sqrt(mu2)) - np.exp(-(x + 1) * np.sqrt(mu2)))
    else:
        F2 = A2 * np.exp(-x * np.sqrt(mu2)) + beta2 / (2 * np.sqrt(mu2) * l2) * \
             (np.exp((1 - x) * np.sqrt(mu1)) - np.exp(-(x + 1) * np.sqrt(mu1)))
    return F2

def _coeffs(T1, S1, K, lambd, theta_or_y, p):
    '''Internal function for Ward and Lough (2011) solution'''
    b11 = T1 * theta_or_y ** 2 + S1 * p + K
    b12 = -K
    b22 = theta_or_y ** 2 + p + K

    mu1 = (b11 / T1 + b22) / 2 + np.sqrt((b11 / T1 + b22) ** 2 / 4 + (b12 ** 2 - b11 * b22) / T1)
    mu2 = (b11 / T1 + b22) / 2 - np.sqrt((b11 / T1 + b22) ** 2 / 4 + (b12 ** 2 - b11 * b22) / T1)
    l1 = T1 + ((mu1 * T1 - b11) / b12) ** 2
    l2 = T1 + ((mu2 * T1 - b11) / b12) ** 2

    beta1 = (mu1 * T1 - b11) / (b12 * 2 * np.pi * p)
    beta2 = (mu2 * T1 - b11) / (b12 * 2 * np.pi * p)

    Delta = 4 * np.sqrt(mu1 * mu2) + 2 * lambd * (np.sqrt(mu1) / l2 + np.sqrt(mu2) / l1)

    A1 = ((lambd / l2 + 2 * np.sqrt(mu2)) * beta1 * np.exp(-np.sqrt(mu1)) - 
           lambd * beta2 / l2 * np.exp(-np.sqrt(mu2))) / Delta / l1

    A2 = (-lambd * beta1 / l1 * np.exp(-np.sqrt(mu1)) + 
          (lambd / l1 + 2 * np.sqrt(mu1)) * beta2 * np.exp(-np.sqrt(mu2))) / Delta / l2

    return b11, b12, b22, mu1, mu2, l1, l2, beta1, beta2, A1, A2

def _safe_factorial(n):
    """Calculate factorial using logarithmic method to avoid overflow."""
    if n < 0:
        return float('inf')
    elif n < 2:
        return 1
    else:
        return np.exp(gammaln(n + 1))

def _StehfestCoeff(jj, N):
    '''Internal function for Ward and Lough (2011) solution'''
    LowerLimit = (jj + 1) // 2
    UpperLimit = min(jj, N // 2)

    V = 0
    for kk in range(LowerLimit, UpperLimit + 1):
        denominator = (_safe_factorial(N // 2 - kk) *
                       _safe_factorial(kk) *
                       _safe_factorial(kk - 1) *
                       _safe_factorial(jj - kk) *
                       _safe_factorial(2 * kk - jj))
        if denominator != 0:  # Prevent division by zero
            V += (kk**(N // 2) * _safe_factorial(2 * kk) / denominator)
    
    V *= (-1) ** (N // 2 + jj)
    return V


# List drawdown and depletion methods so they can be called
# programatically
ALL_DD_METHODS = {'theis': theis,
                    'hunt99ddwn': hunt99ddwn,
                    'wardloughddwn': WardLoughDrawdown}

ALL_DEPL_METHODS = {'glover': glover,
                    'walton': walton,
                    'hunt99': hunt99,
                    'hunt03': hunt2003,
                    'wardlough': WardLoughDepletion}

GPM2CFD = 60*24/7.48 # factor to convert from GPM to CFD