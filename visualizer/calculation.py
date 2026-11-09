# from numba import jit, cuda
import math as m
import numpy as np

# @jit(target_backend='cuda') 
def concentration(drugCon,polymerCon,spaceSize):
    totalSpce=int(spaceSize**3)
    nDrug=int(totalSpce*drugCon)
    monomer=int(totalSpce*polymerCon)
    return nDrug,monomer


# @jit(target_backend='cuda') 
def concentrationPBC(drugCon,polymerCon,spaceSize,polymerSize):
    totalSpce=int(spaceSize**3)
    nDrug=int(totalSpce*drugCon)
    monomer=int(totalSpce*polymerCon)
    nPolymer=int(monomer/polymerSize)
    return nDrug,nPolymer


# @jit(target_backend='cuda') 
def drugPercent(cumDrugRelease,initDrug):
    return (cumDrugRelease/initDrug)*100

# @jit(target_backend='cuda') 
def concentrationProfile(space,spaceSize):
    profile=np.array([])
    for idx in range(spaceSize):
        profile=np.append(profile,np.sum((space[idx,:,:]==2) | (space[idx,:,:]==3)))
    return profile