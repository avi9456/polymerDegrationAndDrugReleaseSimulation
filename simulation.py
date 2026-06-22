
import numpy as np
import randompolymer as polymer
import drug
import calculation as cal
import exportData
import random as rand
def simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,maxPolymerSize,minPolyerSize=1,_saw=False,maxAttempt=100):

    ndrug,nMonomer=cal.concentration(drugCon,polymerCon,spaceSize)

    space = np.zeros((spaceSize,spaceSize,spaceSize),dtype=int)
    data = []
    attempt=maxAttempt
    while(nMonomer and attempt):
        polymerSize=rand.randint(minPolyerSize,maxPolymerSize)
        if polymerSize>=nMonomer:
            polymerSize=nMonomer
        [d,flag]=polymer.randomPolymerGen(polymerSize,spaceSize,space,_saw)
        if not flag:
            attempt-=1
            continue
        data.append(d)
        nMonomer-=polymerSize
        attempt=maxAttempt
    print(len(data))

    drugPos = drug.drugFill(space,ndrug,spaceSize)
    print(np.shape(drugPos))

    MCstep=0
    cumDrugRelease=0
    t_MC=[0,10,100,800]
    snapshots=[{'t_MC':'MCstep','profile':range(0,spaceSize)}]
    MtByMCdata=np.array([['MCStep','M(t)','cumDrugRelease','drugRelease','remainDrug']])
    # print(len(drugPos))
    while MCstep<MC_max and len(drugPos):
        print(MCstep)
        if MCstep in t_MC:
            snapshots.append({'t_MC':MCstep,'profile':cal.concentrationProfile(space,spaceSize)})
        (newPos,drugRelease)=drug.step(drugPos,space,spaceSize,permeablity)
        drugPos=newPos
        remainDrug=len(drugPos)
        cumDrugRelease+=drugRelease
        M_t=cal.M_t(cumDrugRelease,ndrug)
        MtByMCdata=np.append(MtByMCdata,[[MCstep,M_t,cumDrugRelease,drugRelease,remainDrug]],axis=0)
        MCstep+=1
    return MtByMCdata,snapshots
    # exportData.exportData(MtByMCdata,'./data.csv')

if __name__=='__main__':
    spaceSize=25
    drugCon=0.2
    polymerCon=0.6
    permeablity=0.2
    MC_max=1000
    saw=False
    for i in range(50):
        releaseData,snapshots=simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,(spaceSize*3),_saw=saw,maxAttempt=(spaceSize**2))
        exportData.exportData(releaseData,'./report1/releaseData'+str(i)+'.csv')
        snapshots=np.array([[snap['t_MC']]+list(snap['profile']) for snap in snapshots])
        exportData.exportData(snapshots,'./report1/snapshots'+str(i)+'.csv')
