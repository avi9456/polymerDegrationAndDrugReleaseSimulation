# from numba import jit, cuda
import numpy as np
import drug
import simulationFunction as sim
import calculation as cal
import exportData as ed
import matplotlib.pyplot as plt
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

# @jit(target_backend='cuda') 
def worker(para):
    (spaceSize, drugCon, polymerCon, permeablity, MC_max,maxPolymerSize, PBC,maxAttempt,allopen)=para
    releaseData,snapshot,ndrug,trappedDrug=simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,maxPolymerSize,PBC=PBC,maxAttempt=maxAttempt,allOpen=allopen)
    # mt_run=releaseData[1:,0].astype(float)
    remainDrug=releaseData[1:,-1].astype(float)
    cumDrugRelease=releaseData[1:,2].astype(float)
    drugRealease=releaseData[1:,3].astype(float)
    m_t=releaseData[1:,1].astype(float)
    # snapshots.append(snapshot)
    # if len(remainDrug) == 0:
    #     # mt_run= np.zeros(MC_max)
    #     remainDrug = np.full(MC_max, ndrug)
    if len(remainDrug) < MC_max:
        pad_len    = MC_max - len(remainDrug)
        remainDrug = np.pad(remainDrug, (0,pad_len), constant_values=remainDrug[-1])
    if len(cumDrugRelease) < MC_max:
        pad_len    = MC_max - len(cumDrugRelease)
        cumDrugRelease = np.pad(cumDrugRelease, (0,pad_len), constant_values=cumDrugRelease[-1])
    if len(drugRealease) < MC_max:
        pad_len    = MC_max - len(drugRealease)
        drugRealease = np.pad(drugRealease, (0,pad_len), constant_values=0)
    if len(m_t) < MC_max:
        pad_len    = MC_max - len(m_t)
        m_t = np.pad(m_t, (0,pad_len), constant_values=m_t[-1])
    mtStep=np.arange(MC_max)
    releaseData=np.column_stack((mtStep,m_t,cumDrugRelease,drugRealease,remainDrug))
    # mt.append(mt_run)
    # releaseData=np.insert(releaseData,1,remainDrug,axis=1)

    # return (mtStep,m_t,cumDrugRelease,drugRealease,remainDrug)
    return {'releaseData': releaseData, 'ndrug': ndrug, 'trappedDrug': trappedDrug}

def simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,maxPolymerSize,minPolyerSize=2,PBC=False,allOpen=True,maxAttempt=100):
    space = np.zeros((spaceSize,spaceSize,spaceSize),dtype=int)
    if PBC:
        polymerSize=maxPolymerSize
        space,data,drugs,ndrug=sim.fillSpacePBC(space,drugCon,polymerCon,polymerSize,maxAttempt)
    else:
        space,data,drugs,ndrug=sim.fillSpace(space,drugCon,polymerCon,maxPolymerSize,minPolyerSize,maxAttempt)
    releaseData,snapshots,trappedDrug=sim.MC_stepAllOpen(ndrug,drugs,space,spaceSize,permeablity,MC_max,allOpen)
    # print(len(data))
    # print(np.shape(drugPos))
    # ndrug=len(drugPos)

    # MCstep=0
    # cumDrugRelease=0
    # t_MC=[0,10,100,800]
    # snapshots=[{'t_MC':'MCstep','profile':range(0,spaceSize)}]
    # MtByMCdata=np.array([['MCStep','M(t)','cumDrugRelease','drugRelease','remainDrug']])
    # # print(len(drugPos))
    # while MCstep<MC_max and len(drugPos):
    #     print(MCstep)
    #     if MCstep in t_MC:
    #         snapshots.append({'t_MC':MCstep,'profile':cal.concentrationProfile(space,spaceSize)})
    #     (newPos,drugRelease)=drug.step(drugPos,space,spaceSize,permeablity)
    #     drugPos=newPos
    #     remainDrug=len(drugPos)
    #     cumDrugRelease+=drugRelease
    #     M_t=cal.M_t(cumDrugRelease,ndrug)
    #     MtByMCdata=np.append(MtByMCdata,[[MCstep,M_t,cumDrugRelease,drugRelease,remainDrug]],axis=0)
    #     MCstep+=1
    # return MtByMCdata,snapshots
    # exportData.exportData(MtByMCdata,'./data.csv')
    return releaseData,snapshots,ndrug,trappedDrug

if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=4) as executor:
        spaceSize=20
        # drugCon=0.2
        # polymerCon=0.8
        permeablity=0
        MC_max=10000
        PBC=True
        allopen=False
        varyPolymerSize=[5,20,50,100,150,200]
        varyPolymerCon=[0.4,0.5,0.6,0.7]
        varyDrugCon=[0.1,0.15,0.2,0.25,0.3,0.35,0.4]
        dataFolder='./_data/'
        # snapshots=[[]]
        t_axis = np.arange(MC_max)
        count=0
        print('graph1')
        # plt.figure(figsize=(15,10))
        for drugCon in varyDrugCon:
            for polymerCon in varyPolymerCon:
                if (polymerCon+drugCon)>=0.96:
                    break
                # plt.clf()
                psize=0
                
                for maxPolymerSize in varyPolymerSize:
                    psize=maxPolymerSize
                    # mt=[]
                    # drugRemain=[]
                    # for i in range(50):
                    para=[(spaceSize,drugCon,polymerCon,permeablity,MC_max,maxPolymerSize,PBC,(spaceSize**3),allopen) for _ in range(50)]
                        # releaseData,snapshot,ndrug=simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,maxPolymerSize,PBC=PBC,maxAttempt=(spaceSize**3),allOpen=allopen)
                        # mt_run=releaseData[1:,0].astype(float)
                        # remain_run=releaseData[1:,-1].astype(float)
                        # # snapshots.append(snapshot)
                        # if len(mt_run) == 0:
                        #     # mt_run= np.zeros(MC_max)
                        #     remain_run = np.full(MC_max, ndrug)
                        
                        # elif len(mt_run) < MC_max:
                        #     pad_len    = MC_max - len(mt_run)
                        #     # mt_run=np.pad(mt_run,     (0,pad_len), constant_values=mt_run[-1])
                        #     remain_run = np.pad(remain_run, (0,pad_len), constant_values=remain_run[-1])
                        # mt.append(mt_run)
                        # drugRemain.append(remain_run)
                    # averageMt=np.mean(mt,axis=0)
                    # t_axis = np.arange(MC_max)
                    # (drugRemain,releaseData_list,ndrug,trappedDrug,snapshot) = list(executor.map(worker,para))
                    # mtStep,m_t,cumDrugRelease,drugRealease,remainDrug,ndrug,trappedDrug = list(executor.map(worker,para))
                    result = list(executor.map(worker,para))
                    print(np.shape(result[0]['releaseData']))
                    releaseData_list = np.array([res['releaseData'] for res in result])
                    ndrug_list = np.array([res['ndrug'] for res in result])
                    trappedDrug_list = np.array([res['trappedDrug'] for res in result])
                    # trappedDrug_list
                    # for res in result:
                        # np.append(releaseData_list, res['releaseData'])
                        # np.append(ndrug_list, res['ndrug'])
                        # np.append(trappedDrug_list, res['trappedDrug'])
                    print(np.shape(releaseData_list), np.shape(ndrug_list), np.shape(trappedDrug_list))
                    releaseData = np.mean(releaseData_list, axis=0)
                    mean_ndrug = np.mean(ndrug_list)
                    mean_trappedDrug = np.mean(trappedDrug_list)
                    print(np.shape(releaseData),mean_ndrug,mean_trappedDrug)
                    ed.exportDataPandas(releaseData, mean_ndrug, mean_trappedDrug, maxPolymerSize, polymerCon, drugCon, f"{dataFolder}varrypolySize/+exported_data_polymerSize{maxPolymerSize}+polymerCon{polymerCon}+drugCon{drugCon}.json") 
                    # exit()
                    # averageDrugRemain=np.mean(remainDrug,axis=0)

                    # averageDrugRemainpercent=(averageDrugRemain/averageDrugRemain[0])*100
                    # averageDrugRemainpercent = [cal.drugPercent(remain, averageDrugRemain[0]) for remain in averageDrugRemain]
                    # print(f"Average M(t) for polymer size {maxPolymerSize}: {averageMt}")
                    # print(f"Average remaining drug for polymer size {maxPolymerSize}: {averageDrugRemain}")
                
                    # plt.plot(t_axis,averageDrugRemainpercent,label=f'Polymer Size: {maxPolymerSize} | Polymer Con: {polymerCon} | Drug Con: {drugCon}')
        # plt.plot(averageMt,averageDrugRemain)
                # plt.xlabel('MC Steps (time)')
                # plt.ylabel('Remaining Drug (%)')
                # plt.title('Remaining Drug vs Time — Varying Polymer Chain Length')
                # plt.legend()
                # plt.grid(True)
                # plt.tight_layout()
                # plt.savefig(f'{dataFolder}graph1_polymerSize{psize}+polymerCon{polymerCon}+drugCon{drugCon}.png', dpi=150)
                count+=1
                print(count)
        
        plt.clf()
        print('graph2')
        count=0
        # varyPolymerCon=[0.1,0.2,0.3,0.4,0.5,0.6,0.7]
        for drugCon in varyDrugCon:
            for polymerSize in varyPolymerSize:
                # plt.clf() 
                pcon=0  
                for polymerCon in varyPolymerCon:
                    if (polymerCon+drugCon)>=0.96:
                        break
                    pcon=polymerCon
                    # mt=[]
                    # drugRemain=[]
                    # for i in range(50):
                    para=[(spaceSize,drugCon,polymerCon,permeablity,MC_max,polymerSize,PBC,(spaceSize**3),allopen) for _ in range(50)]
                        # releaseData,snapshot,ndrug=simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,polymerSize,PBC=PBC,maxAttempt=(spaceSize**3),allOpen=allopen)
                        # mt_run=releaseData[1:,0].astype(float)
                        # remain_run=releaseData[1:,-1].astype(float)
                        # # snapshots.append(snapshot)
                        # if len(mt_run) == 0:
                        #     mt_run= np.zeros(MC_max)
                        #     remain_run = np.full(MC_max, ndrug)
                        
                        # elif len(mt_run) < MC_max:
                        #     pad_len    = MC_max - len(mt_run)
                        #     # mt_run=np.pad(mt_run,     (0,pad_len), constant_values=mt_run[-1])
                        #     remain_run = np.pad(remain_run, (0,pad_len), constant_values=remain_run[-1])
                        # # mt.append(mt_run)
                        # drugRemain.append(remain_run)
                    # averageMt=np.mean(mt,axis=0)
                    # t_axis = np.arange(MC_max)
                    result = list(executor.map(worker,para))
                    print(np.shape(result[0]['releaseData']))
                    releaseData_list = np.array([res['releaseData'] for res in result])
                    ndrug_list = np.array([res['ndrug'] for res in result])
                    trappedDrug_list = np.array([res['trappedDrug'] for res in result])
                    print(np.shape(releaseData_list), np.shape(ndrug_list), np.shape(trappedDrug_list))
                    releaseData = np.mean(releaseData_list, axis=0)
                    mean_ndrug = np.mean(ndrug_list)
                    mean_trappedDrug = np.mean(trappedDrug_list)
                    print(np.shape(releaseData),mean_ndrug,mean_trappedDrug)
                    ed.exportDataPandas(releaseData, mean_ndrug, mean_trappedDrug, polymerSize, polymerCon, drugCon, f"{dataFolder}varyPolymerCon/exported_data_polymerCon{polymerCon}+polymerSize{polymerSize}+drugCon{drugCon}.json")
                    # drugRemain=list(executor.map(worker,para))
                    # averageDrugRemain=np.mean(drugRemain,axis=0)
                    # averageDrugRemainpercent=(averageDrugRemain/averageDrugRemain[0])*100
                   
                    # averageDrugRemainpercent = [cal.drugPercent(remain, averageDrugRemain[0]) for remain in averageDrugRemain]
                    # print(f"Average M(t) for polymer concentration {polymerCon}: {averageMt}")
                    # print(f"Average remaining drug for polymer concentration {polymerCon}: {averageDrugRemain}")
                
                    # plt.plot(t_axis,averageDrugRemainpercent,label=f'Polymer Concentration: {polymerCon} | drug Concentration: {drugCon} | polymer Size: {polymerSize}')
                # plt.xlabel('MC Steps (time)')
                # plt.ylabel('Remaining Drug (%)')
                # plt.title('Remaining Drug vs Time — Varying Polymer Con Length')
                # plt.legend()
                # plt.grid(True)
                # plt.tight_layout()
                # plt.savefig(f'{dataFolder}graph2_polymerCon{pcon}+polymerSize{polymerSize}+drugCon{drugCon}.png', dpi=150)
                count+=1
                print(count)

        # plt.clf() 
        print('graph3')
        count=0
        for polymerCon in varyPolymerCon:
            for polymerSize in varyPolymerSize:
                # plt.clf()
                dcon=0
                for drugCon in varyDrugCon:
                    if (polymerCon+drugCon)>=0.96:
                        break
                    dcon=drugCon
                    # mt=[]
                    # drugRemain=[]
                    # for i in range(50):
                    para=[(spaceSize,drugCon,polymerCon,permeablity,MC_max,polymerSize,PBC,(spaceSize**3),allopen) for _ in range(50)]
                        # releaseData,snapshot,ndrug=simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,polymerSize,PBC=PBC,maxAttempt=(spaceSize**3),allOpen=allopen)
                        # mt_run=releaseData[1:,0].astype(float)
                        # remain_run=releaseData[1:,-1].astype(float)
                        # # snapshots.append(snapshot)
                        # if len(mt_run) == 0:
                        #     # mt_run= np.zeros(MC_max)
                        #     remain_run = np.full(MC_max, ndrug)
                        
                        # elif len(mt_run) < MC_max:
                        #     pad_len    = MC_max - len(mt_run)
                        #     # mt_run=np.pad(mt_run,     (0,pad_len), constant_values=mt_run[-1])
                        #     remain_run = np.pad(remain_run, (0,pad_len), constant_values=remain_run[-1])
                        # # mt.append(mt_run)
                        # drugRemain.append(remain_run)
                    # averageMt=np.mean(mt,axis=0)
                    # t_axis = np.arange(MC_max)
                    result = list(executor.map(worker,para))
                    print(np.shape(result[0]['releaseData']))
                    releaseData_list = np.array([res['releaseData'] for res in result])
                    ndrug_list = np.array([res['ndrug'] for res in result])
                    trappedDrug_list = np.array([res['trappedDrug'] for res in result])
                    print(np.shape(releaseData_list), np.shape(ndrug_list), np.shape(trappedDrug_list))
                    releaseData = np.mean(releaseData_list, axis=0)
                    mean_ndrug = np.mean(ndrug_list)
                    mean_trappedDrug = np.mean(trappedDrug_list)
                    print(np.shape(releaseData),mean_ndrug,mean_trappedDrug)
                    ed.exportDataPandas(releaseData, mean_ndrug, mean_trappedDrug, polymerSize, polymerCon, drugCon, f"{dataFolder}varydrugCon/exported_data_polymerCon{polymerCon}+polymerSize{polymerSize}+drugCon{drugCon}.json")
                    # drugRemain=list(executor.map(worker,para))
                    # averageDrugRemain=np.mean(drugRemain,axis=0)
                    # averageDrugRemainpercent = [cal.drugPercent(remain, averageDrugRemain[0]) for remain in averageDrugRemain]
                    # averageDrugRemainpercent=(averageDrugRemain/averageDrugRemain[0])*100
                    # print(f"Average M(t) for drug concentration {drugCon}: {averageMt}")
                    # print(f"Average remaining drug for drug concentration {drugCon}: {averageDrugRemain}")
                
                    # plt.plot(t_axis,averageDrugRemainpercent,label=f'Drug Concentration: {drugCon} | Polymer Size: {polymerSize} | Polymer Con: {polymerCon}')
                # plt.xlabel('MC Steps (time)')
                # plt.ylabel('Remaining Drug (%)')
                # plt.title('Remaining Drug vs Time — Varying Drug Concentration')
                # plt.legend()
                # plt.grid(True)
                # plt.tight_layout()
                # plt.savefig(f'{dataFolder}graph3_drugCon{dcon}+polymerSize{polymerSize}+polymerCon{polymerCon}.png', dpi=150)
                count+=1
                print(count)
        # for i in range(50):
        #     releaseData,snapshots=simulation(spaceSize,drugCon,polymerCon,permeablity,MC_max,maxPolymerSize,PBC=PBC,maxAttempt=(spaceSize**2),allOpen=allopen)
        #     exportData.exportData(releaseData,'./report1PBC/releaseData'+str(i)+'.csv')
        #     snapshots=np.array([[snap['t_MC']]+list(snap['profile']) for snap in snapshots])
        #     exportData.exportData(snapshots,'./report1PBC/snapshots'+str(i)+'.csv')
            # print(releaseData)
            # print(snapshots)
