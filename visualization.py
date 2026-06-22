import matplotlib.pyplot as plt
import numpy as np
import randompolymer as polymer
import drug
import calculation as cal
import exportData 
import random as rand
N=15
plt.ion()
ndrug,nMonomerSpace=cal.concentration(0.2,.6,N)
space = np.zeros((N,N,N),dtype=int)
fig = plt.figure(figsize=(20,20))
data = []
# count=10
attempt=100
while(nMonomerSpace and attempt):
    polymerSize=rand.randint(1,N)
    if polymerSize>=nMonomerSpace:
        polymerSize=nMonomerSpace
    [d,flag]=polymer.randomPolymerGen(polymerSize,N,space,True)
    if not flag:
        attempt-=1
        continue
    data.append(d)
    nMonomerSpace-=polymerSize
    attempt=100
print(len(data))
drugPos = drug.drugFill(space,ndrug,N)
print(np.shape(drugPos))
x=[]
y=[]
z=[]
# print(data)
for polydata in data:
    i,j,k = np.hsplit(polydata,3)
    x.append(i)
    y.append(j)
    z.append(k)
# x,y,z=np.hsplit(data,col)
ax=plt.axes(projection='3d')
ax.set_xlim(0, N)
ax.set_ylim(0, N)
ax.set_zlim(0, N)
ax.grid(True)
for i in range(len(x)):
    ax.plot(x[i],y[i],z[i],markersize=50,c='red',alpha=0.5)
    ax.scatter(x[i],y[i],z[i],s=100,c='red')

# print(drugPos)
drugPlot=ax.scatter(drugPos[:,0],drugPos[:,1],drugPos[:,2],s=100,c='blue',alpha=0.7)
# plt.show()
# print(drugPos)
MCstep=150
t_MC=0
cumDrugRelease=0
MtByMCdata=np.array([['MCStep','M(t)','cumDrugRelease','drugRelease','remainDrug']])
# print(len(drugPos))
while MCstep and len(drugPos):
    # plt.show(block=False)
    t_MC+=1
    (newPos,drugRelease)=drug.step(drugPos,space,N,0.2)
    if len(newPos):
        # drugPlot.remove(drugPos[:,0],drugPos[:,1],drugPos[:,2])
        drugPlot.remove()
        drugPlot=ax.scatter(newPos[:,0],newPos[:,1],newPos[:,2],s=100,c='blue',alpha=0.7)
        drugPos=newPos
        remainDrug=len(drugPos)
        # plt.close()
        # plt.show(block=False)
        # plt.draw()
        # print(drugRelease,np.shape(drugPos))

        cumDrugRelease+=drugRelease
        M_t=cal.M_t(cumDrugRelease,ndrug)
        MtByMCdata=np.append(MtByMCdata,[[t_MC,M_t,cumDrugRelease,drugRelease,remainDrug]],axis=0)
    plt.pause(0.01)
    MCstep-=1
exportData.exportData(MtByMCdata,'./data.csv')
plt.ioff()
plt.show()
