import numpy as np
import random as rnd
import directionMap as Map
def drugFill(space,molecul,spaceSize):
    attempt=100
    i=molecul
    drugPos=[]
    while i and attempt:
        x,y,z=rnd.randint(1,spaceSize-2),rnd.randint(1,spaceSize-2),rnd.randint(1,spaceSize-2)
        if space[x,y,z]==0:
            drugPos.append([x,y,z])
            space[x,y,z]=2
            i-=1
            attempt=100
        else:
            attempt-=1
    return np.array(drugPos)


def possibleStep(pos,space,por):
    if space[pos[0],pos[1],pos[2]]==0:
        return (True,False)
    elif space[pos[0],pos[1],pos[2]]==1 and rnd.random()<por:
        return (True,True)
    return False,False


def outOfBound(pos,spaceSize):
    for i in pos:
        if i<0 or i>=spaceSize:
            return True
    return False

# def possibleStep(pos,d,space,por,spaceSize):
#     # x,y,z=Map.map[d]+np.array(pos)
#     x,y,z=pos
#     if space[x,y,z]==0:
#         return (True,pos)
#     x,y,z=Map.map[d]+np.array(pos)
#     if space[x,y,z]==1 and rnd.random()>por:
#         if outOfBound([x,y,z],spaceSize):
#             return (True,[x,y,z])
#         return possibleStep([x,y,z],d,space,por,spaceSize,mov)
#     return (False,0)


# def outOfBound(pos,spaceSize):
#     # if flag:
#     #     newPos=Map.map[d]+np.array(pos)
#     # else:
#     #     newPos=pos
#     for i in pos:
#         if i<0 or i>=spaceSize:
#             return True
#     return False



def step(drugPos,space,spaceSize,por):
    direction=[-1,-2,-3,3,2,1]
    np.random.shuffle(drugPos)
    newPos=[]
    
    drugRelease=0
    for pos in drugPos:
        placed=False
        np.random.shuffle(direction)
        for d in direction:
            x,y,z=Map.map[d] + pos
            if outOfBound([x,y,z],spaceSize):
                drugRelease+=1
                placed=True
                space[pos[0],pos[1],pos[2]]=0
                break
            else: 
                possible,monomer=possibleStep([x,y,z],space,por)
                if possible:
                    if space[pos[0],pos[1],pos[2]]==3:
                        space[pos[0],pos[1],pos[2]]=1
                    else:
                        space[pos[0],pos[1],pos[2]]=0
                    newPos.append([x,y,z])
                    if monomer:
                        space[x,y,z]=3
                    else:
                        space[x,y,z]=2
                    placed=True
                    break
        if not placed:
            newPos.append(pos)
            
    return (np.array(newPos),drugRelease)

            
