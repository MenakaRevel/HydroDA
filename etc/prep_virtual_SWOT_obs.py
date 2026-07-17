#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import math
import errno
import calendar
import datetime
import itertools
import numpy as np
from numpy import ma
from multiprocessing import Pool
import netCDF4 as nc

####################################
# PARAMETERS & PATHS
####################################
def starttime(): return 2000, 1, 1
def endtime(): return 2000, 12, 31
def mapname(): return "conus_06min"
def obs_name(): return "SWOT"
def input_format(): return "nc" # CHANGE THIS to "nc" to read NetCDF files
def obs_dir(): return "/home/revelnil/CaMa-Flood_v4/out/SWOT-conus_06min"
def dam_list(): return "../dat/dam_glb_15min.txt"
def CaMa_dir(): return "/home/revelnil/CaMa-Flood_v4"
def kml_file(): return '../sat/SWOT_Science_sept2015_Nadir.kml'
def mesh_dir(): return "../sat" 

def map_dimension():
    fname = os.path.join(CaMa_dir(), "map", mapname(), "params.txt")
    with open(fname, "r") as f:
        lines = f.readlines()
    nx = int(list(filter(None, lines[0].split(" ")))[0])
    ny = int(list(filter(None, lines[1].split(" ")))[0])
    gsize = float(list(filter(None, lines[3].split(" ")))[0])
    
    # Safely deduce coordinate bounding box (defaults to global if missing)
    west, north = -180.0, 90.0
    if len(lines) > 7:
        try:
            west = float(list(filter(None, lines[4].split(" ")))[0])
            north = float(list(filter(None, lines[7].split(" ")))[0])
        except (IndexError, ValueError):
            pass
            
    return nx, ny, gsize, west, north

def mk_dir(sdir):
    try:
        os.makedirs(sdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise

####################################
# PART 1: 3D MASK PRE-COMPUTATION
####################################
orbit = [] 
target_day = 1 

def calc_dis(latlonlist):
    """Calculates distance for the orbit mesh."""
    global orbit, target_day
    lat, lon = latlonlist[0], latlonlist[1]
    
    diff = np.sqrt((lat - orbit[:, 1])**2 + ((lon - orbit[:, 2]) * math.cos(lat * 3.1416 / 180.))**2)
    min_diff = min(diff[orbit[:, 0] == target_day])
    min_index = int(min(orbit[:, 3][diff == min_diff]))
    
    if min_index == 61068: min_index -= 1
    if min_index == 0: min_index += 1

    x0, y0 = lon, lat
    x1, y1 = orbit[min_index, 2], orbit[min_index, 1]
    x2, y2 = orbit[min_index + 1, 2], orbit[min_index + 1, 1]
    
    if x2 == x1 and y2 == y1:
        x2, y2 = orbit[min_index + 2, 2], orbit[min_index + 2, 1]
        
    x3, y3 = orbit[min_index - 1, 2], orbit[min_index - 1, 1]
    if x3 == x1 and y3 == y1:
        x3, y3 = orbit[min_index - 2, 2], orbit[min_index - 2, 1]

    a1 = (y2 - y1) / ((x2 - x1 + 1e-20) * math.cos(lat * 3.1416 / 180.))
    a2 = (y3 - y1) / ((x3 - x1 + 1e-20) * math.cos(lat * 3.1416 / 180.))

    a = (a1 + a2) / 2.
    b = -1.
    c = 0.
    
    d = abs(a * (x0 - x1) * math.cos(lat * 3.1416 / 180.) + b * (y0 - y1) + c) / math.sqrt(a**2 + b**2)

    ds = -2. / (a1 + a2)
    dsA = d * ds / math.sqrt(1 + ds**2) 
    A = d / math.sqrt(1 + ds**2)

    eqa_len = 40076.954 
    mer_len = 40009.000

    mer_dis = dsA * mer_len / 360.
    eqa_dis = A * eqa_len / 360.

    return math.sqrt(mer_dis**2 + eqa_dis**2)

def dam_loc():
    nx, ny, *_ = map_dimension()
    damloc = np.ones([ny, nx], np.float32)
    with open(dam_list(), "r") as f:
        lines = f.readlines()
    for line in lines[1:]:
        line_parts = list(filter(None, line.split(" ")))
        ix, iy = int(line_parts[4]) - 1, int(line_parts[5]) - 1
        damloc[iy, ix] = -9999.0
    return damloc

def prepare_3d_mask(ncpus):
    """Generates a single 3D array (21, ny, nx) containing final observation masks."""
    global orbit, target_day
    mk_dir(mesh_dir())
    mask_file = os.path.join(mesh_dir(), "swot_mask_21days.bin")
    nx, ny, gsize, west, north = map_dimension()
    
    # Check if unified mask already exists
    if os.path.exists(mask_file):
        print("Unified 21-day mask already exists. Skipping calculation.")
        return

    print("Parsing KML and generating 21-day unified mask...")
    with open(kml_file(), 'r', encoding='latin-1', errors='ignore') as f:
        data = f.readlines()

    incoord, index, day = 0, 0, 1
    for line in data:
        if "Day" in line:
            idx = line.find("Day")
            day = int(line[idx+3 : idx+5])
        if "</coordinates>" in line: incoord = 0
        if incoord == 1:
            parts = line.split(",")
            orbit.append([day, float(parts[1]), float(parts[0]), index])
            index += 1
        if "<coordinates>" in line: incoord = 1

    orbit = np.asarray(orbit)
    
    # DYNAMICALLY generate coordinates based EXACTLY on your map size
    latlist = np.arange(north - gsize/2, north - (ny * gsize), -gsize)[:ny]
    lonlist = np.arange(west + gsize/2, west + (nx * gsize), gsize)[:nx]
    latlonlist = list(itertools.product(latlist, lonlist))
    
    # Load static features for masking
    rivwdth_thr = 50.0 
    wth_file = os.path.join(CaMa_dir(), "map", mapname(), "rivwth_gwdlr.bin")
    rivwth = np.fromfile(wth_file, np.float32).reshape(ny, nx)
    dams = 1.0 # dam_loc()
    
    # Pre-allocate 3D mask array
    mask_3d = np.zeros((21, ny, nx), dtype=np.float32)
    p = Pool(ncpus)

    for day_iter in range(1, 22):
        target_day = day_iter
        print(f"Calculating and masking Day {target_day}/21...")
        
        outmesh = p.map(calc_dis, latlonlist)
        mesh_in = np.asarray(outmesh, dtype=np.float32).reshape([ny, nx])
        
        # Determine SWOT swath (10km to 60km off-nadir)
        SWOTmesh = (mesh_in >= 10) & (mesh_in <= 60)
        
        # Create final binary mask logic (must be within swath, river width, and dam checks)
        final_mask = (SWOTmesh & (rivwth >= rivwdth_thr) & (dams > 0.0)).astype(np.float32)
        mask_3d[day_iter - 1, :, :] = final_mask
        
    p.terminate()
    
    # Save the unified array
    print(f"Saving unified 3D mask to {mask_file}...")
    mask_3d.tofile(mask_file)

####################################
# PART 2: OBSERVATION ERROR
####################################
def calc_static_obs_error():
    nx, ny, *_ = map_dimension()
    k, q, ovs_err = 1.00, 1.00, 0.10
    
    len_file = os.path.join(CaMa_dir(), "map", mapname(), "rivlen.bin")
    wth_file = os.path.join(CaMa_dir(), "map", mapname(), "rivwth_gwdlr.bin")
    next_file = os.path.join(CaMa_dir(), "map", mapname(), "nextxy.bin")
    
    rivlen = np.fromfile(len_file, np.float32).reshape(ny, nx) * 0.0 + 1.0 
    rivwth = np.fromfile(wth_file, np.float32).reshape(ny, nx) * 1.0e-3
    nextx = (np.fromfile(next_file, np.int32).reshape(2, ny, nx)[0] != -9999) * 1.0
    
    area = (k * rivlen) * (q * rivwth)
    
    obs_err = ovs_err * (1 / (k * rivlen + 1.0e-20)) * (1 / (q * rivwth + 1.0e-20)) * nextx
    obs_err = obs_err * (area >= 1.0) * 1.0 + 0.25 * (1 / (k * rivlen + 1.0e-20)) * (1 / (q * rivwth + 1.0e-20)) * nextx * (area < 1.0) * 1.0
    obs_err = ma.masked_where(area < 0.625, obs_err).filled(0.25) 
    obs_err = obs_err * ((obs_err <= 0.25) * 1.0) + 0.25 * ((obs_err > 0.25) * 1.0)
    obs_err = obs_err * nextx
    return obs_err.astype(np.float32)

def get_days(syear, eyear):
    return (datetime.date(eyear, 12, 31) - datetime.date(syear, 1, 1)).days + 1

####################################
# PART 3: VIRTUAL OBSERVATION PIPELINE
####################################
def prepare_single_obs_set(indir, outdir):
    mk_dir(outdir)
    nx, ny, *_ = map_dimension()
    syear, smon, sday = starttime()
    eyear, emon, eday = endtime()
    file_fmt = input_format()
    
    print("Pre-calculating static observation error field...")
    static_obs_err = calc_static_obs_error()
    
    print("Loading 21-day unified mask into RAM...")
    mask_file = os.path.join(mesh_dir(), "swot_mask_21days.bin")
    ALL_21_MASKS = np.fromfile(mask_file, np.float32).reshape(21, ny, nx)
    
    total_days = get_days(syear, eyear)
    SWOTOBS = np.zeros([total_days, ny, nx], np.float32)
    
    indays = 0
    start_dt = datetime.date(syear, 1, 1)
    
    print(f"Processing surface elevation ({file_fmt} format) and applying noise...")
    for year in range(syear, eyear + 1):
        nt = 366 if calendar.isleap(year) else 365
        
        # Format parsing logic
        if file_fmt == "bin":
            fname = os.path.join(indir, f"sfcelv{year}.bin")
            yr_data = np.fromfile(fname, np.float32).reshape([nt, ny, nx])
        elif file_fmt == "nc":
            fname = os.path.join(indir, f"o_sfcelv{year}.nc")
            with nc.Dataset(fname, 'r') as ds:
                yr_data = ds.variables['sfcelv'][:].astype(np.float32)
        else:
            raise ValueError("Invalid input_format. Choose 'bin' or 'nc'.")
            
        SWOTOBS[indays:indays+nt, :, :] = yr_data
        
        for day in range(nt):
            target_dt = datetime.date(year, 1, 1) + datetime.timedelta(days=day)
            
            # Apply Noise
            np.random.seed(indays + day)
            noise_array = np.random.normal(0.0, static_obs_err, (ny, nx))
            SWOTOBS[indays + day, :, :] += noise_array
            
            # Instantly grab the mask from RAM instead of disk
            days_diff = (target_dt - datetime.date(syear, smon, sday)).days
            swot_day_index = days_diff % 21
            obs_mask = ALL_21_MASKS[swot_day_index, :, :]
            
            SWOTOBS[indays + day, :, :] = ma.masked_where(obs_mask != 1.0, SWOTOBS[indays + day, :, :]).filled(-9999.0)
            
        indays += nt
        
    print("Computing mean and standard deviation over time...")
    masked_obs = ma.masked_equal(SWOTOBS, -9999.0)
    SWOTmean = np.mean(masked_obs, axis=0)
    SWOTstd = np.std(masked_obs, axis=0)
    
    print("Writing TXT output files...")
    for day in range(total_days):
        target_dt = start_dt + datetime.timedelta(days=day)
        yyyy, mm, dd = f"{target_dt.year:04d}", f"{target_dt.month:02d}", f"{target_dt.day:02d}"
        
        txtfile = os.path.join(outdir, f"sfcelv_{yyyy}{mm}{dd}.txt")
        with open(txtfile, "w") as txtf:
            for ix in range(nx):
                for iy in range(ny):
                    val = SWOTOBS[day, iy, ix]
                    if val != -9999.0:
                        txtf.write(f"{ix:04d}  {iy:04d}    {val:10.4f}  {SWOTmean[iy, ix]:10.4f}  {SWOTstd[iy, ix]:10.4f}  SWOT\n")

####################################
# MAIN EXECUTION
####################################
if __name__ == "__main__":
    ncpus = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    
    # 1. Generate/Ensure the single 3D mask file exists
    prepare_3d_mask(ncpus)
    
    # 2. Process exactly ONE setup
    indir = "/home/revelnil/CaMa-Flood_v4/out/SWOT-conus_06min"
    outdir = "../obs"
    
    print(f"Starting pipeline from {indir} -> {outdir}")
    prepare_single_obs_set(indir, outdir)
    print("Done.")