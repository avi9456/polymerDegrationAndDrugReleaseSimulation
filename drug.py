# from numba import jit, cuda
import numpy as np
import random as rnd
import directionMap as Map


class drugMolecule:
    def __init__(self,postion:np.array,traped:bool=True) -> None:
        self.postion=np.array(postion)
        self.traveled=set((tuple(postion)))
        self.traped=traped
        self.moves=0
    def isTraped(self,moves)->bool:
        if len(self.traveled)>moves and self.traped:
            self.traped=False
        else:
            self.traveled.clear()
        return self.traped

    def possibleStep(self,targetPos,space,por)-> bool:
        if space[targetPos[0],targetPos[1],targetPos[2]]==0:
            return True
        elif space[targetPos[0],targetPos[1],targetPos[2]]==1 and rnd.random()<por:
            return True
        return False
    
    def outOfBound(self,targetPos,spaceSize:int)-> bool:
        for i in targetPos:
            if i<0 or i>=spaceSize:
                return True
        return False
    def step(self,space,por,allOpen=True)->bool:
        direction=[-1,-2,-3,3,2,1]
        spaceSize=len(space)
        np.random.shuffle(direction)
        for d in direction:
            x,y,z=Map.map[d] + np.array(self.postion)
            if not allOpen:
                y,z=(np.array([y,z])%spaceSize)
            if self.outOfBound([x,y,z],spaceSize):
                space[self.postion[0],self.postion[1],self.postion[2]]=0
                self.traped=False
                return True
            else: 
                possible=self.possibleStep([x,y,z],space,por)
                if possible:
                    if space[self.postion[0],self.postion[1],self.postion[2]]==3:
                        space[self.postion[0],self.postion[1],self.postion[2]]=1
                    else:
                        space[self.postion[0],self.postion[1],self.postion[2]]=0
                    self.postion=[x,y,z]
                    if space[x,y,z]==1:
                        space[x,y,z]=3
                    else:
                        space[x,y,z]=2
                    if self.traped:
                        self.traveled.add(tuple(self.postion))
                    return False
        return False

def drugFill(space,nDrug,maxAttempt=100):
    attempt=maxAttempt
    drugs=[]
    while nDrug and attempt:
        if attempt<=0:
            nDrug-=1
            attempt=maxAttempt
        x,y,z=rnd.randint(1,len(space)-2),rnd.randint(1,len(space)-2),rnd.randint(1,len(space)-2)
        if space[x,y,z]==0:
            space[x,y,z]=2
            nDrug-=1
            attempt=maxAttempt
            drugs.append(drugMolecule(np.array([x,y,z]),True))
        else:
            attempt-=1
    return drugs

# @jit(target_backend='cuda') 
# old drugFill function
# def drugFill(space,molecul,spaceSize,maxAttempt=100):
#     attempt=maxAttempt
#     i=molecul
#     drugPos=[]
#     while i and attempt:
#         x,y,z=rnd.randint(1,spaceSize-2),rnd.randint(1,spaceSize-2),rnd.randint(1,spaceSize-2)
#         if space[x,y,z]==0:
#             drugPos.append([x,y,z])
#             space[x,y,z]=2
#             i-=1
#             maxAttempt-=1
#             attempt=maxAttempt
#         else:
#             attempt-=1
#     return np.array(drugPos)

# old possibleStep function
# # @jit(target_backend='cuda') 
# def possibleStep(pos,space,por):
#     if space[pos[0],pos[1],pos[2]]==0:
#         return (True,False)
#     elif space[pos[0],pos[1],pos[2]]==1 and rnd.random()<por:
#         return (True,True)
#     return False,False

# # @jit(target_backend='cuda') 
# old outOfBound function
# def outOfBound(pos,spaceSize):
#     for i in pos:
#         if i<0 or i>=spaceSize:
#             return True
#     return False

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


# @jit(target_backend='cuda') 
# old step function
# def step(drugPos,space,spaceSize,por,allOpen=True):
#     direction=[-1,-2,-3,3,2,1]
#     np.random.shuffle(drugPos)
#     newPos=[]
    
#     drugRelease=0
#     for pos in drugPos:
#         placed=False
#         np.random.shuffle(direction)
#         for d in direction:
#             x,y,z=Map.map[d] + pos
#             if not allOpen:
#                 y,z=(np.array([y,z])%spaceSize)
#             if outOfBound([x,y,z],spaceSize):
#                 drugRelease+=1
#                 placed=True
#                 space[pos[0],pos[1],pos[2]]=0
#                 break
#             else: 
#                 possible,monomer=possibleStep([x,y,z],space,por)
#                 if possible:
#                     if space[pos[0],pos[1],pos[2]]==3:
#                         space[pos[0],pos[1],pos[2]]=1
#                     else:
#                         space[pos[0],pos[1],pos[2]]=0
#                     newPos.append([x,y,z])
#                     if monomer:
#                         space[x,y,z]=3
#                     else:
#                         space[x,y,z]=2
#                     placed=True
#                     break
#         if not placed:
#             newPos.append(pos)
            
#     return (np.array(newPos),drugRelease)

            
