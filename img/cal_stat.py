import numpy as np
from numpy import ma
import math
import scipy as sp

import params as pm
#--
def rmse(assimilations,truths,add=1.0e-20):
  rms = []
  assimilations = np.array(assimilations)
  truths        = np.array(truths)
  # assimilations is 2d arry
  for i in range(0,len(truths)): 
    rms.append(np.sqrt(((assimilations[:,i]-truths[i]*np.ones([pm.ens_mem()],np.float32))**2).mean()))
  return rms#(rms/abs(truths+add))*100.0 
#--
def AI(assimilations,corruptions,truths):
  AIval = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)
  #--
  error =abs(truths-np.mean(corruptions,axis=1))/(truths+1e-20) 
  AIval = 1. - abs((truths-np.mean(assimilations,axis=1))/(truths-np.mean(corruptions,axis=1))+1.0e-20)

  return ma.masked_where(error<0.1,AIval).filled(0.0)
#--
def AI_new(assimilations,corruptions,truths):
  AIval = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)
  #--
  error =abs(truths-np.mean(corruptions,axis=1))/(truths+1e-20)  
  AIval = 1. - abs((np.mean(assimilations,axis=1)-np.mean(corruptions,axis=1))/(truths-np.mean(corruptions,axis=1)+1.0e-20)-1.)
 
  #return ma.masked_where(error<0.1,AIval).filled(0.0) 
  return AIval,error
#--
def AI_new1(assimilations,corruptions,truths):
  AIval = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)
  #--
  error =abs(truths-np.mean(corruptions,axis=1))/(truths+1e-20)
#  AIval =  abs((truths-np.mean(assimilations,axis=1))/(truths-np.mean(corruptions,axis=1)+1.0e-20)-1.)
  AIval =abs(1. - (np.mean(assimilations,axis=1)-np.mean(corruptions,axis=1))/(truths-np.mean(corruptions,axis=1)+1.0e-20))

  return ma.masked_where(error<0.1,AIval).filled(0.0)
#--
def pBias(assimilations,corruptions,truths):
  pB = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)

  #pB   = 100.0 * np.nansum(np.nanmean(assimilations,axis=1)-truths)/(np.nansum(truths)+1e-20)
  #pB_c = 100.0 * np.nansum(np.nanmean(corruptions,axis=1)-truths)/(np.nansum(truths)+1e-20)

  pB   = 100.0 * np.sum(np.mean(assimilations,axis=1)-truths)/(np.sum(truths)+1e-20)
  pB_c = 100.0 * np.sum(np.mean(corruptions,axis=1)-truths)/(np.sum(truths)+1e-20)

  #pB   = 100.0 * (np.nanmean(assimilations,axis=1)-truths/(np.nanmean(truths)+1e-20))
  #pB_c = 100.0 * (np.nanmean(corruptions,axis=1)-truths/(np.nanmean(truths)+1e-20))

  return pB,pB_c
#-- Refernce Oubanas_etal2018
def RMSE(assim,true):
    """Root Mean Squre Error"""
    asm=np.mean(assim,axis=1)
    org=np.array(true)
    T=float(len(org))
    RootMSE=math.sqrt((1/T)*sum(asm-org)**2)
    return RootMSE
#---
def rRMSE(assim,true):
    """relative Root Mean Square Error"""
    asm=np.mean(assim,axis=1)
    org=np.array(true)
    T=float(len(org))
    rRootMSE=math.sqrt((1/T)*sum((asm-org)/(org+1.0e-20))**2)
    return rRootMSE
#---
def NRMSE(assim,true):
    """Normalized Root Mean Square Error"""
    asm=np.mean(assim,axis=1)
    org=np.array(true)
    orgQ=np.mean(org)
    T=float(len(org))
    NRootMSE=(1/(orgQ+1.0e-20))*math.sqrt((1/T)*sum(asm-org)**2)
    return NRootMSE

#----
def VE(assim,true):
    """Volumetric Efficency"""
    asm=np.mean(assim,axis=1)
    org=np.array(true)
    orgQ=np.mean(org)
    T=float(len(org))
    VolEff=1.0 - (sum(np.abs(asm-org))/sum(org))
    return VolEff
#----
def NSE(assims,corrs,true):
  assims = np.array(assims)
  assim  = np.mean(assims,axis=1)
  corrs  = np.array(corrs)
  corr   = np.mean(corrs,axis=1)
  true   = np.array(true)
  #--
  Qm = np.mean(true)
  #---
  NSA=1.0 - (np.sum((assim-true)**2)/(sum((Qm-true)**2)+1.0e-20))
  NSC=1.0 - (np.sum((corr-true)**2)/(sum((Qm-true)**2)+1.0e-20))
  #--
  NSEc=((NSA-NSC)/(1.0-NSC))

  return NSEc, NSA, NSC
#----
def PRI(assims,corrs,true,tol=15):
  assims = np.array(assims)
  assim  = np.mean(assims,axis=1)
  corrs  = np.array(corrs)
  corr   = np.mean(corrs,axis=1)
  true   = np.array(true)
  #--
  Qt_max = np.amax(true)
  Ttp    = np.where(true==Qt_max)[0][0]
  Tp0    = max(Ttp-tol,0)
  Tp1    = min(Ttp+tol,len(true))
  #---
  Qa_max = np.amax(assim[Tp0:Tp1])
  Qc_max = np.amax(corr[Tp0:Tp1])
  #---
  Ttp    = float(Ttp)
  Tap    = float(np.where(assim==Qa_max)[0][0])
  Tcp    = float(np.where(corr==Qc_max)[0][0])
  #--
  PDRI   = 1.0 -(abs(Qt_max-Qa_max)/(abs(Qt_max-Qc_max)+1.0e-20))
  PTRI   = 1.0 -(abs(Ttp-Tap)/(abs(Ttp-Tcp)+1.0e-20))
  #--
  return PDRI, PTRI
#----
def EnsSpr(assim,corr,true):
    assim = np.array(assim)
    corr  = np.array(corr)
    true  = np.array(true)
    #---
    assim_max= np.amax(assim,axis=1)
    assim_min= np.amin(assim,axis=1)
    corr_max = np.amax(corr ,axis=1)
    corr_min = np.amin(corr ,axis=1)
    #---
    ENSPR=1.0 - ((assim_max-assim_min)/(corr_max-corr_min))
    ENSPR=ma.masked_less(ENSPR,0.0).filled(0.0)
    return ENSPR
#----
def AQ(PDRI,PTRI,meanAI,meanEnsSpr,w1=0.3,w2=0.3,w3=0.4,w4=0.2):
    return w1*PDRI + w2*PTRI + w3*meanAI + w4*meanEnsSpr
#----
def KGE(assims,corrs,true):
  """Kling etal 2009, Kling et al 2011"""
  assims = np.array(assims)
  assim  = np.mean(assims,axis=1)
  corrs  = np.array(corrs)
  corr   = np.mean(corrs,axis=1)
  true   = np.array(true)
  #----
  CC=sp.stats.pearsonr(assim,true)
  BR=np.mean(assim)/np.mean(true)
  RV=(np.mean(assim)/np.std(assim))/(np.mean(true)/np.std(true))

  return 1 - math.sqrt((CC-1)**2+(BR-1)**2+(RV-1)**2)
