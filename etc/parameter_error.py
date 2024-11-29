#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import scipy
import scipy.signal
from scipy.stats import truncnorm, lognorm
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
def spatially_correlated_random(nx,ny,mu=1.0,sigma=1.0,tau=50.0):
    """
    Spatially correlated random variable 
    spital dimension are [nx,ny]
    tau_s = decorrelation length
    """
    correlation_scale = tau
    x = np.arange(-correlation_scale, correlation_scale)
    y = np.arange(-correlation_scale, correlation_scale)
    X, Y = np.meshgrid(x, y)
    dist = np.sqrt(X*X + Y*Y)
    filter_kernel = np.exp(-dist**2/(2*correlation_scale))
    # noise = np.random.lognormal(mu, sigma, size=(ny,nx)) #.reshape(ny,nx)
    noise = np.random.normal(mu, sigma, size=(ny,nx))
    noise = scipy.signal.fftconvolve(noise, filter_kernel, mode='same')
    # noise = np.real(noise)*(1.0/(np.std(np.real(noise))+1e-20))
    noise = np.real(noise)*(1.0/(np.max(np.abs(noise))+1e-20))
    # noise = np.fft.fft2(noise)
    # filter_kernel = np.fft.fft2(filter_kernel)
    # noise = np.fft.ifft2(noise*filter_kernel)
    # noise = np.real(noise)*(1.0/(np.std(np.real(noise))+1e-20))
    return noise
#====================================================================================
def truncated_lognormal(mu,sigma,threshold,N=1):
    """
    Truncated lognormal distribution
    """
    
    # Calculate the parameters of the truncated distribution
    alpha = (threshold - mu) / sigma
    truncation_point = (threshold - mu) / sigma

    # Create a truncated lognormal distribution
    trunc_lognorm = truncnorm(a=0.0, b=threshold, loc=mu, scale=sigma)

    # Generate random samples from the truncated lognormal distribution
    errors = trunc_lognorm.rvs(size=N)

    return errors
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

    # errors = errors/(np.std(errors)+1e-20)

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

    # Convert to float32
    flp=np.float32(flp)

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

    # Error of the river depth - regional bias (ß)
    beta = spatially_correlated_random(nx,ny,mu=1.0,sigma=1.12,tau=1000.0)*1.12 # std of 1.12

    # Random error (ε)
    # efs = err_wd(0.0,1.12,N=nx*ny).reshape(ny,nx)#*1.12 # std of 1.12
    efs = truncated_lognormal(0.0,1.12,2.5,N=nx*ny).reshape(ny,nx)#*1.12 # std of 1.12

    # Add the error to the river depth
    dph_ = dph * beta + efs #+ err_wd(1.0,1.12,N=1) # 39% error
    
    # Remove extremely large values and low values
    dph_[dph_>np.max(dph)] = np.max(dph)
    dph_[dph_<1.0] = 1.0

    # Convert to float32
    dph_=np.float32(dph_)

    # Save the corrupted floodplain map
    dph_.tofile(CaMa_dir+"/map/"+mapname+"/rivhgt_corrupt.bin")
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
    wth = np.fromfile(CaMa_dir+"/map/"+mapname+"/rivwth.bin",np.float32).reshape(ny,nx)

    # Error of the river depth - regional bias (ß)
    beta = spatially_correlated_random(nx,ny,mu=1.0,sigma=1.22,tau=1000.0)*1.22 # std of 1.12

    # Random error (ε)
    # efs = err_wd(1.0,1.0,N=nx*ny).reshape(ny,nx)*1.22 # std of 1.12 
    efs = truncated_lognormal(0.0,1.22,5.0,N=nx*ny).reshape(ny,nx)#*1.12 # std of 1.12

    # Add the error to the river depth
    wth_ = wth * beta + efs #+ err_wd(1.0,1.22,N=1) # 39% error

    # Remove extremely large values and low values
    wth_[wth_>np.max(wth)] = np.max(wth)
    wth_[wth_<5.0] = 5.0

    # Convert to float32
    wth_=np.float32(wth_)

    # Save the corrupted floodplain map
    wth_.tofile(CaMa_dir+"/map/"+mapname+"/rivwth_corrupt.bin")
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
    man = np.fromfile(CaMa_dir+"/map/"+mapname+"/rivman.bin",np.float32).reshape(ny,nx)

    # Random error (ε)
    efs = err_wd(0.0,1.69,N=nx*ny).reshape(ny,nx) #*1.69 # std of 1.12 

    # Add the error to the river depth
    man_ = man * efs #+ err_n(1.0,1.69,N=1) # 39% error

    # Remove extremely large values and low values
    man_[man_>0.035] = 0.035
    man_[man_<0.025] = 0.025

    # Convert to float32
    man_=np.float32(man_)

    # Save the corrupted floodplain map
    man_.tofile(CaMa_dir+"/map/"+mapname+"/rivman_corrupt.bin")
    print ("save ----> "+CaMa_dir+"/map/"+mapname+"/rivman_corrupt.bin")
    return 0
#====================================================================================
if __name__=="__main__":
    # corrupt_flddph("glb_15min") # Corrupt the global 15min resolution floodplain map
    # corrupt_rivhgt("glb_15min") # Corrupt the global 15min resolution river depth map
    corrupt_rivwth("glb_15min") # Corrupt the global 15min resolution river width map
    # corrupt_rivman("glb_15min") # Corrupt the global 15min resolution river Manning's n map