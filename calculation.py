import math as m
import numpy as np


def concentration(drugCon,polymerCon,spaceSize):
    totalSpce=int(spaceSize**3)
    nDrug=int(totalSpce*drugCon)
    monomer=int(totalSpce*polymerCon)
    return nDrug,monomer


def M_t(cumDrugRelease,initDrug):
    return (cumDrugRelease/initDrug)*100


def concentrationProfile(space,spaceSize):
    profile=np.array([])
    for idx in range(spaceSize):
        profile=np.append(profile,np.sum((space[idx,:,:]==2) | (space[idx,:,:]==3)))
    return profile