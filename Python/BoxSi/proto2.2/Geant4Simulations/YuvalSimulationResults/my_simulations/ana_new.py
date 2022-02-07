# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 16:24:51 2021

@author: Yuval Was
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib qt
N=4




def merge_pulses(indices,times,window,energies,thresh):
    pulses=dict()
    times=np.array(list(map(float,times)))
    energies=np.array(list(map(float,energies)))
    indices=np.array(list(map(int,indices)))
    
    
    for i in range(N):
        pulses[i]=np.column_stack((times[indices==i], energies[indices==i]))
        pulses[i]=pulses[i][pulses[i][:,0].argsort()]

    
    for key in pulses.keys():
        temp_time=np.array([])
        temp_energy=np.array([])
        det_pulses=pulses[key]
        for i in range(det_pulses.shape[0]):
            if i==0:
                temp_time=np.append(temp_time,det_pulses[0,0])
                temp_energy=np.append(temp_energy,det_pulses[0,1])
            elif (det_pulses[i,0]-temp_time[-1])>window:
                temp_time=np.append(temp_time,det_pulses[i,0])
                temp_energy=np.append(temp_energy,det_pulses[i,1])
            else:
                temp_energy[-1]+=det_pulses[i,1]
        t_e_mat=np.column_stack((temp_time, temp_energy))
        pulses[key]=t_e_mat[t_e_mat[:,1]>thresh]
        
        """
        print(pulses[key])
        print("-----")
    print("-@@-@@-@@-@@-@@-")"""
    return pulses



c_names=["Time Gamma","Energy Gamma", "Id Gamma","Time Neutron","Energy Neutron","Id Neutron"]
df = pd.read_csv('rossi_nt_rossi.csv',sep=',',   header=9, names=c_names)

energy_threshold=[0.4,2.6] #Energy Threshold for Gamma detection and for Neutron Detection in MeV
rate=595500 #rate of fission events per second

t_global=0;
particles_names=["Gamma","Neutron"]
det_id=[]
final_pulses_times=np.array([])
final_pulses_real_times=np.array([])
final_pulses_energy=np.array([])
final_pulses_id=np.array([])
final_pulses_isN=np.array([])

for i_row in range(df.shape[0]):
    if (i_row%10000==0):
        print(i_row)
    t_global+=np.random.exponential(1/rate); #FIX HERE, account for unmeasure events.
    for particle_i in range(0,2):
        particle_name=particles_names[particle_i]
        row=df.iloc[i_row]
        energy=row["Energy "+particle_name]
        
        if not isinstance(energy,str):
            continue
        if not isinstance(row["Time "+particle_name],str):
            continue
        energy=np.array(list(map(float,row["Energy "+particle_name].split(";"))))
        det_index=np.array(list(map(float,row["Id "+particle_name].split(";"))))
        time=np.array(list(map(float,row["Time "+particle_name].split(";"))))
    
        index_keep = np.array([energy[i]>energy_threshold[particle_i] for i in range(len(energy))])
        energy=energy[index_keep]
        det_index=det_index[index_keep]
        time = time[index_keep]
        
        final_pulses_times=np.concatenate([final_pulses_times,time+t_global])
        final_pulses_real_times=np.concatenate([final_pulses_real_times,time])
        final_pulses_energy=np.concatenate([final_pulses_energy,energy])
        final_pulses_id=np.concatenate([final_pulses_id,det_index])
        final_pulses_isN=np.concatenate([final_pulses_isN,particle_i*np.ones(time.shape)])

final_pulses_times=final_pulses_times*1e9 #to ns

titles=["Time (ns)","Energy(MeV)","Detector ID", "IsNeutron"]
df_final=pd.DataFrame(np.column_stack((final_pulses_times,final_pulses_energy,final_pulses_id,final_pulses_isN)),columns=titles)
df_final=df_final.sort_values(by="Time (ns)")
df_final.to_csv("pulse_train.csv")

plt.figure()
plt.hist(final_pulses_energy)
plt.figure()
plt.hist(final_pulses_id)
plt.show()


            
            
    
    
    