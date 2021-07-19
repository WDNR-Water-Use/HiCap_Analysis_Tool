import numpy as np
import pandas as pd
from pathlib import Path
import pytest


datapath = Path('hicap_analysis/tests/data')

@pytest.fixture
def theis_results():
    excel_file = datapath / 'HighCap_Analysis_Worksheet_Example.xlsm'
    params = pd.read_excel(excel_file,
        sheet_name='Property_Drawdown_Analysis',
        usecols='C:D',
        skiprows=7,
        index_col=0)

    theis_res = pd.read_excel(excel_file, 
        sheet_name='Theis', 
        skiprows= 9,
        usecols=('A:B,H:I'),
        names=['well1_dd','well1_r','well2_dd','well2_r'])
    return {'params':params, 'theis_res':theis_res}

def test_theis(theis_results):
    from hicap_analysis import wells as wo
    # read the two Q values and convert to CFD
    Q = [float(i)*wo.GPM2CFD for i in [theis_results['params'].loc['Pumping Rate Well #1 (gpm)'],
            theis_results['params'].loc['Pumping Rate Well #2 (gpm)']]]
    S = float(theis_results['params'].loc['Storage Coefficient (unitless)'])
    T = float(theis_results['params'].loc['Transmissivity (ft2/day)'])
    #dist = [float(i) for i in [theis_results['params'].loc['Distance from Well #1 to Municpal Well'],
    #        theis_results['params'].loc['Distance from Well #2 to Municpal Well']]]
    dist = theis_results['theis_res'].well1_r
    time = float(theis_results['params'].loc['Theis Time of Pumping (days)'])
    
    dd = [wo._theis(T, S, time, dist, currQ) for currQ in Q]
    assert np.allclose(dd[0],theis_results['theis_res'].well1_dd, atol=0.5)
    assert np.allclose(dd[1],theis_results['theis_res'].well2_dd, atol=0.5)