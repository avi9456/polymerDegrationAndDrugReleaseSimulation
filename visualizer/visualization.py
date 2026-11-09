# from numba import jit, cuda
import matplotlib.pyplot as plt
import numpy as np
import polymer as polymer
import drug
import calculation as cal 
import random as rand
import simulationFunction as sim

# @jit(target_backend='cuda') 
def visualization(spaceSize,drugCon,polymerCon,permeability,MC_max,maxPolymerSize,minPolyerSize=1,fixPolymer=False,maxAttempt=1000):
    plt.ion()
    space = np.zeros((spaceSize,spaceSize,spaceSize),dtype=int)
    fig = plt.figure(figsize=(10,10))
    if fixPolymer:
        polymerSize=maxPolymerSize
        space,data,drugPos,ndrug=sim.fillSpacePBC(space,drugCon,polymerCon,polymerSize,maxAttempt)
    else:
        space,data,drugPos,ndrug=sim.fillSpace(space,drugCon,polymerCon,maxPolymerSize,minPolyerSize,maxAttempt)
    # ndrug,nMonomerSpace=cal.concentration(drugCon,polymerCon,spaceSize)
    # ndrug=len(drugPos)
    # data = []
    # attempt=maxAttempt
    # while(nMonomerSpace and attempt):
    #     polymerSize=rand.randint(minPolyerSize,maxPolymerSize)
    #     if polymerSize>=nMonomerSpace:
    #         polymerSize=nMonomerSpace
    #     [d,flag]=polymer.randomPolymerGen(polymerSize,spaceSize,space,_saw)
    #     if not flag:
    #         attempt-=1
    #         continue
    #     data.append(d)
    #     nMonomerSpace-=polymerSize
    #     attempt=maxAttempt
    # print(len(data))
    # drugPos = drug.drugFill(space,ndrug,spaceSize)
    # print(np.shape(drugPos))
    x=[]
    y=[]
    z=[]
    # print(data)
    for polydata in data:
        i,j,k = np.hsplit(polydata,3)
        x.append(i)
        y.append(j)
        z.append(k)

    ax=plt.axes(projection='3d')
    ax.set_xlim(0, spaceSize)
    ax.set_ylim(0, spaceSize)
    ax.set_zlim(0, spaceSize)
    ax.grid(True)
    for i in range(len(x)):
        ax.scatter(x[i],y[i],z[i],s=100,c='red')
        ax.plot(x[i],y[i],z[i],markersize=50,c='red',alpha=0.5)
        

    drugPlot=ax.scatter(drugPos[:,0],drugPos[:,1],drugPos[:,2],s=100,c='blue',alpha=0.7)
    
    MCstep=0
    # print(len(drugPos))
    while MCstep < MC_max and len(drugPos):
        # plt.show(block=False)
        # t_MC+=1
        (newPos,drugRelease)=drug.step(drugPos,space,spaceSize,permeability)
        if len(newPos):
            drugPlot.remove()
            drugPlot=ax.scatter(newPos[:,0],newPos[:,1],newPos[:,2],s=100,c='blue',alpha=0.7)
            drugPos=newPos
        plt.pause(0.01)
        MCstep+=1


if __name__=='__main__':
    spaceSize=15
    drugCon=0.2
    polymerCon=0.6
    permeablity=0.0
    MC_max=1000
    saw=True
    visualization(spaceSize,drugCon,polymerCon,permeablity,MC_max,(spaceSize*3),maxAttempt=(spaceSize**2),fixPolymer=True)
