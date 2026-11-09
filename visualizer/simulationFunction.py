# from numba import jit, cuda
import numpy as np
import polymer as polymer
import drug
import calculation as cal
import random as rand

# @jit(target_backend='cuda') 
def fillSpace(space,drugCon,polymerCon,maxPolymerSize,minPolyerSize=1,maxAttempt=100):
    spaceSize=space.shape[0]
    ndrug,nMonomer=cal.concentration(drugCon,polymerCon,spaceSize)
    data = []
    attempt=maxAttempt
    while(nMonomer and attempt):
        polymerSize=rand.randint(minPolyerSize,maxPolymerSize)
        if polymerSize>=nMonomer:
            polymerSize=nMonomer
        [d,flag]=polymer.randomPolymerGen(polymerSize,spaceSize,space)
        if not flag:
            attempt-=1
            continue
        data.append(d)
        nMonomer-=polymerSize
        maxAttempt-=polymerSize
        attempt=maxAttempt
    # print(len(data))
    drugPos = drug.drugFill(space,ndrug,spaceSize,maxAttempt)
    # print(np.shape(drugPos))
    return space,data,drugPos,ndrug

# @jit(target_backend='cuda') 
def fillSpacePBC(space,drugCon,polymerCon,polymerSize=15,maxAttempt=100):
    spaceSize=space.shape[0]
    ndrug,nPolymer=cal.concentrationPBC(drugCon,polymerCon,spaceSize,polymerSize)
    data = []
    attempt=maxAttempt
    while(nPolymer and attempt):
        [d,flag]=polymer.randomPolymerGenPBC(polymerSize,spaceSize,space,)
        if not flag:
            attempt-=1
            continue
        data.append(d)
        nPolymer-=1
        maxAttempt-=polymerSize
        attempt=maxAttempt
    # print(len(data))
    drugPos = drug.drugFill(space,ndrug,spaceSize,maxAttempt)
    print(np.shape(drugPos))
    return space,data,drugPos,ndrug
# @jit(target_backend='cuda') 
def MC_stepAllOpen(ndrug,drugPos,space,spaceSize,permeablity,MC_max,allOpen=True):

    MCstep=1
    cumDrugRelease=0
    t_MC=[0,10,100,800]
    snapshots=[{'t_MC':'MCstep','profile':range(0,spaceSize)}]
    MtByMCdata=np.array([['MCStep','drugReleasePercent','cumDrugRelease','drugRelease','remainDrug']])
    # print(len(drugPos))
    while MCstep<=MC_max and len(drugPos):
        # print(MCstep)
        if MCstep in t_MC:
            snapshots.append({'t_MC':MCstep,'profile':cal.concentrationProfile(space,spaceSize)})
        (newPos,drugRelease)=drug.step(drugPos,space,spaceSize,permeablity,allOpen)
        drugPos=newPos
        remainDrug=len(drugPos)
        cumDrugRelease+=drugRelease
        M_t=cal.drugPercent(cumDrugRelease,ndrug)
        MtByMCdata=np.append(MtByMCdata,[[MCstep,M_t,cumDrugRelease,drugRelease,remainDrug]],axis=0)
        MCstep+=1
    return MtByMCdata,snapshots