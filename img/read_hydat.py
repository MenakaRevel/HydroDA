# usr/lib/python 

import numpy as np
import re
import shutil
import os
import datetime
import sys
import pandas as pd
import sqlite3
#--
def get_hydat_loc(rivername,fname='/work/a06/menaka/Hydat/data/Hydata_cmf_list.csv'):
    # read Hydata_cmf_list.csv
    df=pd.read_csv(fname)
    # print (df)
    df=df.loc[df['River'].str.contains(rivername),:]
    
    # print (df['Hydat_id'].values,df['Station'].values,df['ix1'].values,df['iy1'].values)
    #
    return (list( df['Hydat_id'].values),
    list(df['Station'].values),
    list(df['ix1'].values),
    list(df['iy1'].values))
#--
def get_hdat_station(station,fname='/work/a06/menaka/Hydat/data/Hydata_cmf_list.csv'):
    # read Hydata_cmf_list.csv
    df=pd.read_csv(fname)
    # print (df)
    df=df[df['Station'].str.contains(station)]
    # print (df)
    #
    return df['ix1'].values[0],df['iy1'].values[0],df['ix2'].values[0],df['iy2'].values[0]
#--
def hydat_dis(station,syear,eyear,smon=1,emon=12,sday=1,eday=31):
    # initlize dataframe
    Date = pd.date_range(start=datetime.date(syear,smon,sday).strftime('%Y-%m-%d'), 
    end=datetime.date(eyear,emon,eday).strftime('%Y-%m-%d'), freq="D")
    flowdata = pd.DataFrame(
        np.full((len(Date), 2), -9999.0), columns=["FLOW", "QC"], index=Date
    )
    CA_HYDAT='/work/a06/menaka/Hydat/data/Hydat.sqlite3'
    con = sqlite3.connect(CA_HYDAT)
    sqlstat = "select * from DLY_FLOWS WHERE STATION_NUMBER = ?"
    Readed_Streamflow = pd.read_sql_query(sqlstat, con, params=[station])
    Readed_Streamflow = Readed_Streamflow[Readed_Streamflow["YEAR"] >= syear]
    Readed_Streamflow = Readed_Streamflow[Readed_Streamflow["YEAR"] <= eyear]

    if len(Readed_Streamflow) != 0:
            ### loop read streamflow data
        for index, row in Readed_Streamflow.iterrows():
            NDays = row["NO_DAYS"]
            for iday in range(1, NDays + 1):
                cdate = pd.to_datetime(
                    {"year": [row["YEAR"]], "month": [row["MONTH"]], "day": [iday]}
                ).values
                #            cdates = pd.to_datetime(str(row['YEAR'])+'-'+str(row['MONTH'])+'-'+str(iday))
                if (
                    row["FLOW" + str(iday)] != np.nan
                    and row["FLOW" + str(iday)] != None
                    and float(row["FLOW" + str(iday)]) > 0
                ):
                    flowdata.loc[cdate, "FLOW"] = row["FLOW" + str(iday)]
                    flowdata.loc[cdate, "QC"] = row["FLOW" + "_SYMBOL" + str(iday)]

    return flowdata['FLOW'].values