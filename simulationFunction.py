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
    drugs = drug.drugFill(space,ndrug,maxAttempt)
    # print(np.shape(drugPos))
    return space,data,drugs,ndrug

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
    drugs = drug.drugFill(space,ndrug,maxAttempt)
    
    # print(np.shape(drugs))
    return space,data,drugs,ndrug
# @jit(target_backend='cuda') 
def MC_stepAllOpen(ndrug,drugs,space,spaceSize,permeablity,MC_max,allOpen=True):

    MCstep=1
    cumDrugRelease=0
    t_MC=[0,10,100,800]
    snapshots=[{'t_MC':'MCstep','profile':range(0,spaceSize)}]
    MtByMCdata=np.array([['MCStep','drugReleasePercent','cumDrugRelease','drugRelease','remainDrug']])
    # print(len(drugPos))
    while MCstep<=MC_max and len(drugs):
        # print(MCstep)
        tmpDrugs=[]
        if MCstep in t_MC:
            snapshots.append({'t_MC':MCstep,'profile':cal.concentrationProfile(space,spaceSize)})
        for drug in drugs:
            if not drug.step(space,permeablity,allOpen):
                if((MCstep+1)%100==0):
                    drug.isTraped(30)
                tmpDrugs.append(drug)
            else:
                del drug
        # (newPos,drugRelease)=drug.step(drugs,space,spaceSize,permeablity,allOpen)
        # drugPos=newPos
        # remainDrug=len(drugPos)
        drugRealease = len(drugs)-len(tmpDrugs)
        drugs=tmpDrugs
        remainDrug=len(drugs)
        cumDrugRelease+=drugRealease
        M_t=cal.drugPercent(cumDrugRelease,ndrug)
        MtByMCdata=np.append(MtByMCdata,[[MCstep,M_t,cumDrugRelease,drugRealease,remainDrug]],axis=0)
        MCstep+=1
    trappedDrug=0
    for drug in drugs:
        if drug.isTraped(30):
            trappedDrug+=1
    return MtByMCdata,snapshots,trappedDrug