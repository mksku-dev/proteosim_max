from pyteomics import achrom

def predict_lc_retention_times(peptides):
    """
    Predicts the retention time (rounded to two decimals) of each peptide and retruns a dict of each peptide with its retetnion time.

    Parameters
    ----------
    peptides: list 
        List of peptide sequence. 

    Returns
    -------
    dict 
        Dictonary of each peptide with its retention time.
    """
    result = {}
    for peptide in peptides:
        rt = achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0)
        result[peptide] = round(rt, 2)
    return result  

import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Plots the retention time as a histogramm. Meaning it counts the amount of proteins per retetention time.

    Parameters
    ----------
    resolution : int
        Gives the amount of bins inside the histogramm (default = 30)
    
    retention_times: 

    Returns
    -------
    plot : 
        Shows a Histogramm with the retention times a y- scale and the peptides as x-scale
    """
    plt.hist(list(retention_times), bins=resolution)

    plt.xlabel('Retention time (min)')
    plt.ylabel('Number of peptides')

    plt.show()


def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Looks at the retention time dictionary and filters out every peptide which is not in the time window
    Parameters
    ----------
    peptide_rt_map: dict
        Dictonary of every Peptide with its corresponding Retention time value 
    lower_ret_time :  float
        Minimum peptide retention time to keep
    upper_ret_time : float
        Maximum peptide retention time to keep

    Returns
    -------
    list of str
        Returns every peptide which retention time is in the time window

    """
   
    filtered_peptides = [pep for pep in peptide_rt_map if lower_ret_time <= peptide_rt_map[pep] <= upper_ret_time]
    return filtered_peptides
    # dont have the retention times but should have !!! ig but also dont it doesnt matter 
