import calendar
import pandas as pd
#from hicap_analysis.wells import GPM2CFD
import numpy as np

def Q2ts(depl_pump_time, depletion_years, Q_scalar):
    """Function to convert a scalar Q to a time series

    Args:
        depl_pump_time (float): Number of days the well should be pumping from day 0
        depletion_years (float): Number of years of pumping to be simulated
        Q_scalar (float): Pumping rate in GPM

    Returns:
        Q (pd.Series): day-indexed time series of pumping in CFD
    """
    y1 = np.zeros(365)
    y1[:depl_pump_time] = Q_scalar
    Q = pd.Series(index = range(1,(depletion_years*365)+1),
                data = list(y1)*depletion_years) * GPM2CFD
    return Q

def create_timeseries_template(filename='../examples/blank_ts.csv', 
                numyears=1, 
                well_ids = ['well1','well2']):
    """create a template timeseries CSV file ready to populate with 
    pumping rates in gallons per minute for running multiple times

    Args:
        filename (str, optional): Filename in which to write the template. Defaults to '../examples/blank_ts.csv'.
        numyears (int, optional): Number of repeated years to generate. Defaults to 1.
        well_ids(list, optional): list of well identifiers which will be used as column names
    """
    
    # do not set to true, unless maybe you are running a single year
    leapyear = False
    if leapyear:
        refyear = 2020
    else:
        refyear = 2021
        
    # make lists of months and days for the date columns
    months = []
    days = []
    for i in range(1,13):
        dayspermonth = calendar.monthrange(refyear,i)[1]
        months += ([i]*dayspermonth)
        days += list(range(1,dayspermonth+1))

    # now make a dataframe
    ts_df = pd.DataFrame(data = {'year':sum([[i]*365 for i in range(1,numyears+1)],[]),
                                            'month':months * numyears,
                                            'day':days*numyears})
    ts_df.index += 1
    ts_df['sequential_day'] = ts_df.index
    for cc in well_ids:
        ts_df[cc] = 0
    # and write it out
    ts_df.to_csv(filename, index=None)    
if __name__=="__main__":
    create_timeseries_template(well_ids=[f'well{i}' for i in range(1,5)])