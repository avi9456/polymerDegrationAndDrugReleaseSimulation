# from numba import jit, cuda
import numpy as np
import random
import directionMap as Map

# @jit(target_backend='cuda') 
def insidesize(postion,size):
    # print(postion)
    for p in postion:
        if p<0 or p>=size:
            return False
    return True

# @jit(target_backend='cuda') 
def cleanSpace(polymer,space):
    for cord in polymer:
        space[cord[0],cord[1],cord[2]]=0

# @jit(target_backend='cuda') 
def randomPolymerGen(chainLen,size,space):
    x,y,z = random.randint(0,size-1),random.randint(0,size-1),random.randint(0,size-1)
    if space[x,y,z]>0:
        return [np.array([]),False]
    polymer=np.array([[x,y,z]])
    space[x,y,z]=1
    # print (polymer)
    chainLen-=1
    while(chainLen):
        possible_d=[-3,-2,-1,1,2,3]
        random.shuffle(possible_d)
        monomerPlaced=False
        for d in possible_d:
            x,y,z=Map.map[d]+polymer[-1]
            if insidesize([x,y,z],size) and space[x,y,z]==0:
                # print('inside size')
                
                polymer=np.append(polymer,[[x,y,z]],axis=0)
                space[x,y,z]=1
                chainLen-=1
                monomerPlaced=True
                break
        if not monomerPlaced:
            cleanSpace(polymer,space)
            return [np.array([]),False]
    # if not saw:
    #     for cord in polymer:
    #         space[cord[0],cord[1],cord[2]]=1
    return (polymer,True)


# @jit(target_backend='cuda') 
def randomPolymerGenPBC(chainLen,size,space):
    x,y,z = random.randint(0,size-1),random.randint(0,size-1),random.randint(0,size-1)
    if space[x,y,z]>0:
        return [np.array([]),False]
    polymer=np.array([[x,y,z]])
    space[x,y,z]=1
    chainLen-=1
    monomerPlaced=False
    while(chainLen):
        possible_d=[-3,-2,-1,1,2,3]
        random.shuffle(possible_d)
        monomerPlaced=False
        for d in possible_d:
            x,y,z=(Map.map[d]+polymer[-1])%size
            if x<0:
                x+=size
            elif y<0:
                y+=size
            elif z<0:
                z+=size
            if space[x,y,z]==0:
                polymer=np.append(polymer,[[x,y,z]],axis=0)
                chainLen-=1
                monomerPlaced=True
                space[x,y,z]=1
                break
        if not monomerPlaced:
            cleanSpace(polymer,space)
            return (np.array([]),False)
    return (polymer,True) 

# @jit(target_backend='cuda') 
def randomWalk(polymerChain,space,spaceSize,permeability):
    newPos=[]
    monomerPos=0
    for cord in polymerChain:
        possible_d=[-3,-2,-1,1,2,3]
        random.shuffle(possible_d)
        monomerMoved=False
        for d in possible_d:
            x,y,z=Map.map[d]+cord
            if insidesize([x,y,z],spaceSize):
                if space[x,y,z]==0:
                    newPos.append(np.array([x,y,z]))
                    monomerMoved=True
                    break
                if random.random()<permeability:
                    newPos.append(np.array([x,y,z]))
                    monomerMoved=True
                    break



# print(randomPolymerGen(12,7))