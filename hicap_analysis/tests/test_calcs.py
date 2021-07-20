from os import pardir
import numpy as np
import pandas as pd
from pathlib import Path
import pytest


datapath = Path('hicap_analysis/tests/data')

@pytest.fixture
def theis_results():
    from hicap_analysis import wells as wo
    excel_file = datapath / 'HighCap_Analysis_Worksheet_Example.xlsm'
    p = pd.read_excel(excel_file,
        sheet_name='Property_Drawdown_Analysis',
        usecols='C:D',
        skiprows=7,
        index_col=0)
    # read the two Q values and convert to CFD
    Q = [float(i)*wo.GPM2CFD for i in [p.loc['Pumping Rate Well #1 (gpm)'],
            p.loc['Pumping Rate Well #2 (gpm)']]]
    S = float(p.loc['Storage Coefficient (unitless)'])
    T = float(p.loc['Transmissivity (ft2/day)']) 
    time = float(p.loc['Theis Time of Pumping (days)'])
    params = {'Q':Q,'S':S,'T':T, 'time':time}
    theis_res = pd.read_excel(excel_file, 
        sheet_name='Theis', 
        skiprows= 9,
        usecols=('A:B,H:I'),
        names=['well1_dd','well1_r','well2_dd','well2_r'])
    
    walton_res = pd.read_excel(excel_file,
        sheet_name='Stream#1_Depletion', 
        skiprows= 104,
        usecols=('C,M:N,R,AB:AC,AK'),
        names=['twell','dep1','dep2','timage','rch1','rch2', 'total_dep'])
      
    return {'params':params, 'theis_res':theis_res}


@pytest.fixture
def walton_results():
    from hicap_analysis import wells as wo
    excel_file = datapath / 'HighCap_Analysis_Worksheet_Example.xlsm'  
    walton_res = pd.read_excel(excel_file,
        sheet_name='Stream#1_Depletion', 
        skiprows= 104,
        usecols=('C,M:N,R,AB:AC,AK'),
        names=['t_well','dep1','dep2','t_image','rch1','rch2', 'total_dep'])
    p = pd.read_excel(excel_file,
        sheet_name='Stream#1_Depletion', 
        skiprows= 70,
        nrows=30,
        usecols=('B:D'),
        names=['par','v1','v2'],
        index_col=0,)
    Q = p.loc['Q - Pumping rate (ft^3/dy)'].T.values
    S = p.loc['S - Storage (unitless)'].T.values
    dist = p.loc['a - Distance (feet)'].T.values
    T_gpd_ft = p.loc['T - Transmissivity (gpd/ft)'].T.values
    # a little trickery to get the index of the time array for start and end of each well
    Q_start_day = [pd.to_datetime(i).day_of_year - 1 
        for i in p.loc['  First Day of Annual Pumping ='].T.values]
    Q_end_day = [pd.to_datetime(i).day_of_year - 1 
        for i in p.loc['  Last Day of Annual Pumping ='].T.values]
    params = {'Q':Q, 'S':S, 'T_gpd_ft':T_gpd_ft,
            'Q_start_day':Q_start_day,
            'Q_end_day':Q_end_day, 'dist':dist}
  
  
    return {'params':params, 'walton_res':walton_res}

def test_theis(theis_results):
    """Test for the theis calculations - compared with two wells at multiple distances
        in the example spreadsheet

    Args:
        theis_results (@fixture, dict): parameters and results from example spreadsheet
    """


    from hicap_analysis import wells as wo
    pars = theis_results['params']
    dist = theis_results['theis_res'].well1_r
    
    time = pars['time']
    dd = [wo._theis(pars['T'], pars['S'], time, dist, currQ) for currQ in pars['Q']]
    assert np.allclose(dd[0],theis_results['theis_res'].well1_dd, atol=0.5)
    assert np.allclose(dd[1],theis_results['theis_res'].well2_dd, atol=0.5)

def test_glover(theis_results):
    """Athens test for the glover calculations

    Args:
        theis_results (@fixture, dict): parameters and results from example spreadsheet
    """
    from hicap_analysis import wells as wo
    pars = theis_results['params']
    dist = theis_results['theis_res'].well1_r
    Qs = wo._glover(pars['T'], pars['S'], pars['time'], dist, pars['Q'][0])
    assert all(np.isnan(Qs)== False)
   
def test_sdf():
    """Test for streamflow depletion factor
        using values from original Jenkins (1968) paper
        https://doi.org/10.1111/j.1745-6584.1968.tb01641.x
        note Jenkins rounded to nearest 10 (page 42)
    """
    from hicap_analysis import wells as wo
    dist = 5280./2.
    T = 5.0e4/7.48
    S = 0.5
    sdf = wo._sdf(T,S,dist)
    assert np.allclose(sdf, 520, atol=1.5)

def test_walton(walton_results):
    """Test of a single year to be sure the Walton calculations are made correctly

    Args:
        walton_results ([type]): [description]
    """
    from hicap_analysis import wells as wo

    res = walton_results['walton_res']
    pars = walton_results['params']
    
    dep={}
    rch={}
    for idx in [0,1]:
        dep[idx] = wo._walton(pars['T_gpd_ft'][idx],
                    pars['S'][idx],
                    pars['dist'][idx],
                    res.t_well,
                    pars['Q'][idx]
                    )
        rch[idx] = wo._walton(pars['T_gpd_ft'][idx],
                    pars['S'][idx],
                    pars['dist'][idx],
                    res.t_image,
                    pars['Q'][idx]
                    )
    dep_tot = dep[0]-rch[0] + dep[1]-rch[1]
    assert np.allclose(dep[0], res.dep1)
    assert np.allclose(dep[1], res.dep2)
    assert np.allclose(rch[0], -res.rch1)
    assert np.allclose(rch[1], -res.rch2)
    assert np.allclose(dep_tot, res.total_dep)