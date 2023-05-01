#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import re
"""
Prepare errors for CaMa-Flood parameters.
"""
#====================================================================================
#-----------------------------------------------------------------------------------
# Parameter   | Error Equation    |           Errors         |          Reference
#------------------------------------------------------------------------------------
# Floodplain  | flp'= flp + ßflp | Normal error of 2m [MERIT]| Yamazaki et al. (2017)
# Width       | w' = w x ßw x εw | logN(1,1.12^2)  or 39%    | Andreadis et al. (2013)
# Depth       | d' = d x ßd x εd | logN(1,1.22^2)  or 33%    | Moody and Troutman (2002)
# Manning's n | n' = n x εn      | logN(1,1.69^2)  or 50%    | Anees et al. (2017)
#====================================================================================
# functions for error modelling
#====================================================================================
def err_flp(bias_mean,bias_std,elen=10):
    """
    Error for floodplain.
    Calculate the perturbed hypsometric curve by adding the regional bias to the elevation counts.
    bias_mean - Mean of the Gaussian distribution
    bias_std  - Standard deviation of the Gaussian distribution
    elen      - Number of elevation classes
    """
    # Define the regional bias as a Gaussian normal distribution
    regional_bias = np.random.normal(loc=bias_mean, scale=bias_std, size=elen)

    return regional_bias
#====================================================================================
def err_wd(mu,sigma,N=1):
    """
    Error for width and depth.
    Calculate the perturbed lognormally distributed errors.
    mu     - Mean of the lognormal distribution
    sigma  - Standard deviation of the lognormal distribution
    N      - Number of random numbers to generate
    """
    # Generate the lognormally distributed errors
    errors = np.random.lognormal(mu, sigma, N)

    return errors
#====================================================================================
def err_n(mu, sigma, N=1):
    """
    Error for Mannings's n.
    Calculate the perturbed lognormally distributed errors.
    mu     - Mean of the lognormal distribution
    sigma  - Standard deviation of the lognormal distribution
    N      - Number of random numbers to generate
    """
    # Generate the lognormally distributed errors
    errors = np.random.lognormal(mu, sigma, N)

    return errors
#====================================================================================
def corrupt_flddph(mapname,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"):
    """
    Corrupt the floodplain parameters of CaMa-Flood.
    """
    fname=CaMa_dir+"/map/"+mapname+"/params.txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    #-------
    nx     = int(filter(None, re.split(" ",lines[0]))[0])
    ny     = int(filter(None, re.split(" ",lines[1]))[0])
    gsize  = float(filter(None, re.split(" ",lines[3]))[0])

    # Load the original floodplain map
    flp = np.fromfile(CaMa_dir+"/map/"+mapname+"/fldhgt.bin",np.float32).reshape(10,ny,nx)

    # Add the regional bias to the elevation counts
    # flp = flp + err_flp(0,2,elen=10) # 2m error
    err = err_flp(0, 2, elen=10)  # 2m error
    flp = flp + np.broadcast_to(err[:,np.newaxis, np.newaxis], (flp.shape[0], ny, nx))

    # keep consistency with the original floodplain map
    flp[1:] = np.maximum(flp[1:], flp[:-1])
    flp[0]  = np.maximum(flp[0], 0.0)

    # Save the corrupted floodplain map
    flp.tofile(CaMa_dir+"/map/"+mapname+"/fldhgt_corrupt.bin")
    print ("save ----> "+CaMa_dir+"/map/"+mapname+"/fldhgt_corrupt.bin")
    return 0
#====================================================================================
def corrupt_rivhgt(mapname,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"):
    """
    Corrupt the river depth parameters of CaMa-Flood.
    """
    fname=CaMa_dir+"/map/"+mapname+"/params.txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    #-------
    nx     = int(filter(None, re.split(" ",lines[0]))[0])
    ny     = int(filter(None, re.split(" ",lines[1]))[0])
    gsize  = float(filter(None, re.split(" ",lines[3]))[0])

    # Load the original floodplain map
    dph = np.fromfile(CaMa_dir+"/map/"+mapname+"/rivhgt.bin",np.float32).reshape(ny,nx)

    # Add the error to the river depth
    dph = dph + err_wd(1.0,1.12,N=1) # 39% error

    # Save the corrupted floodplain map
    dph.tofile(CaMa_dir+"/map/"+mapname+"/rivhgt_corrupt.bin")
    print ("save ----> "+CaMa_dir+"/map/"+mapname+"/rivhgt_corrupt.bin")
    return 0
#====================================================================================
def corrupt_rivwth(mapname,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"):
    """
    Corrupt the river width parameters of CaMa-Flood.
    """
    fname=CaMa_dir+"/map/"+mapname+"/params.txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    #-------
    nx     = int(filter(None, re.split(" ",lines[0]))[0])
    ny     = int(filter(None, re.split(" ",lines[1]))[0])
    gsize  = float(filter(None, re.split(" ",lines[3]))[0])

    # Load the original floodplain map
    dph = np.fromfile(CaMa_dir+"/map/"+mapname+"/rivwth.bin",np.float32).reshape(ny,nx)

    # Add the error to the river depth
    dph = dph + err_wd(1.0,1.22,N=1) # 39% error

    # Save the corrupted floodplain map
    dph.tofile(CaMa_dir+"/map/"+mapname+"/rivwth_corrupt.bin")
    print ("save ----> "+CaMa_dir+"/map/"+mapname+"/rivwth_corrupt.bin")
    return 0
#====================================================================================
def corrupt_rivman(mapname,CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"):
    """
    Corrupt the river width parameters of CaMa-Flood.
    """
    fname=CaMa_dir+"/map/"+mapname+"/params.txt"
    with open(fname,"r") as f:
        lines=f.readlines()
    #-------
    nx     = int(filter(None, re.split(" ",lines[0]))[0])
    ny     = int(filter(None, re.split(" ",lines[1]))[0])
    gsize  = float(filter(None, re.split(" ",lines[3]))[0])

    # Load the original floodplain map
    dph = np.fromfile(CaMa_dir+"/map/"+mapname+"/rivman.bin",np.float32).reshape(ny,nx)

    # Add the error to the river depth
    dph = dph + err_n(1.0,1.69,N=1) # 39% error

    # Save the corrupted floodplain map
    dph.tofile(CaMa_dir+"/map/"+mapname+"/rivman_corrupt.bin")
    print ("save ----> "+CaMa_dir+"/map/"+mapname+"/rivman_corrupt.bin")
    return 0
#====================================================================================
if __name__=="__main__":
    # corrupt_flddph("glb_15min") # Corrupt the global 15min resolution floodplain map
    # corrupt_rivhgt("glb_15min") # Corrupt the global 15min resolution river depth map
    # corrupt_rivwth("glb_15min") # Corrupt the global 15min resolution river width map
    corrupt_rivman("glb_15min") # Corrupt the global 15min resolution river Manning's n map