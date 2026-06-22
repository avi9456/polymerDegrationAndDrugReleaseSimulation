import numpy as np
import random
import directionMap as Map

def insidesize(postion,size):
    # print(postion)
    for p in postion:
        if p<0 or p>=size:
            return False
    return True
def cleanSpace(polymer,space):
    for cord in polymer:
        space[cord[0],cord[1],cord[2]]=0
def randomPolymerGen(chainLen,size,space,saw=False):
    x,y,z = random.randint(0,size-1),random.randint(0,size-1),random.randint(0,size-1)
    if space[x,y,z]>0:
        return [np.array([]),False]
    polymer=np.array([[x,y,z]])
    if saw:
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
                if saw:
                    space[x,y,z]=1
                chainLen-=1
                monomerPlaced=True
                break
        if not monomerPlaced:
            if saw:
                cleanSpace(polymer,space)
            return [np.array([]),False]
    if not saw:
        for cord in polymer:
            space[cord[0],cord[1],cord[2]]=1
    return (polymer,True)



    


# print(randomPolymerGen(12,7))