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
    return {'params':params, 'theis_res':theis_res}

def test_theis(theis_results):
    """Test for the theis calculations - compared with two wells at multiple distances
        in the example spreadsheet

    Args:
        theis_results ([@fixture, dict]): parameters and results from example spreadsheet
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
        theis_results ([@fixture, dict]): parameters and results from example spreadsheet
    """
    pars = theis_results['params']
    dist = theis_results['theis_res'].well1_r
    from hicap_analysis import wells as wo
    Qs = wo._glover(pars['T'], pars['S'], pars['time'], dist, pars['Q'][0])
    assert all(np.isnan(Qs)== False)
    j=2