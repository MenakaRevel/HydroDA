#!/usr/bin/python
import os
import re
import datetime
import numpy as np 
import json
import pandas as pd
import requests
from io import StringIO
################################
# Functions to read SWOT data
################################
def get_swot_loc(rivername,mapname="Mackenzie_06min",fname="/cluster/data6/menaka/HydroDA/dat/SWOT_alloc_Mackenzie_06min.txt"):
    # read metadata
    meta="/work/a06/menaka/SWOT/SWOT_Mackenzie_Station_list.txt"
    dfmeta=pd.read_csv(meta,sep=",")
    names=dfmeta[dfmeta["basin"]==rivername]["identifier"].values
    #read fname
    dfcgls=pd.read_csv(fname,delim_whitespace=True)
    pname=dfcgls[dfcgls["station"].isin(names)]["station"].values
    xlist=dfcgls[dfcgls["station"].isin(names)]["ix"].values - 1
    ylist=dfcgls[dfcgls["station"].isin(names)]["iy"].values - 1
    egm08=dfcgls[dfcgls["station"].isin(names)]["EGM08"].values
    egm96=dfcgls[dfcgls["station"].isin(names)]["EGM96"].values
    return pname,xlist,ylist,egm08,egm96
################################
def query_hydrocron(station, start_time="2024-01-01", end_time="2024-11-01"):
    """Query Hydrocron for reach-level time series data.

    Parameters
    ----------
    station: str - String SWORD reach identifier
    start_time: str - String time to start query format: "2024-01-01T00:00:00Z"
    end_time: str - String time to end query

    Returns
    -------
    pandas.DataFrame that contains query results
    """
    query_url= "https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries"
    fields = "node_id,time_str,wse,wse_u,wse_r_u,node_q"
    params = {
        "feature": 'Node',
        "feature_id": str(station),
        "output": "csv",
        "start_time": start_time+"Z",
        "end_time": end_time+"Z",
        "fields": fields
    }
    results = requests.get(query_url, params=params)
    # print (results, params)
    if "results" in results.json().keys():
        results_csv = results.json()["results"]["csv"]
        df = pd.read_csv(StringIO(results_csv))
        # Remove fill values for missing observations
        df = df.loc[(df["wse"] != -999999999999.0)].reset_index(drop=True)
        df = df.loc[(df['node_q']<=1) & (df['wse_r_u']<1.0), :]
        # Convert time_str to datetime format
        df.time_str = pd.to_datetime(df.time_str)
        # do a outlier removal
        # Remove outliers from 'wse' column using 1.5 * IQR
        Q1, Q3 = df["wse"].quantile([0.25, 0.75])
        IQR = Q3 - Q1
        df = df[(df["wse"] >= Q1 - 1.5 * IQR) & (df["wse"] <= Q3 + 1.5 * IQR)]
        df = df.loc[:,['time_str','wse']]
    else:
        df = pd.DataFrame({
            "time_str": [datetime.datetime(1900, 1, 1).strftime("%Y-%m-%dT%H:%M:%S")],
            "wse": [np.nan],
            })

    return df
################################
def swot_WSE(station,syear,eyear,smon=1,emon=10,sday=1,eday=31):
    start = datetime.datetime(syear, smon, sday,0,0,0)
    end = datetime.datetime(eyear, emon, eday,23,59,59)
    # time = (end - start).days + 1  # time in days
    # data = np.full(time, -9999.0, dtype=np.float32)  # WSE in [m], -9999.0 for no observations

    df=query_hydrocron(station,start.strftime("%Y-%m-%dT%H:%M:%S"),end.strftime("%Y-%m-%dT%H:%M:%S"))

    df["time_str"] = pd.to_datetime(df["time_str"], format="%Y-%m-%dT%H:%M:%S")

    # Filter and update data efficiently using DataFrame operations
    # mask = (df["datetime"] >= start) & (df["datetime"] <= end)
    # df_filtered = df[mask].copy()
    df["days_from_start"] = (df["time_str"] - start).dt.days

    # Separate data and time into arrays
    time_array = df["days_from_start"].values
    data_array = df["wse"].values #+ egm08 - egm96

    return time_array, data_array