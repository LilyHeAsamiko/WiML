# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
pd.options.display.float_format = '${:,.8f}'.format
from sklearn import preprocessing
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import unicodedata
import os
import scipy.special as sps
from scipy.special import erf

path = r'D:\TUT\Medical\biophysics\experiment\single cell sequencing'
#importing the dataset
#dataset2 = np.loadtxt('D:/TUT/Medical/biophysics/experiment/TF predictivity.txt', delimiter ='\t',dtype = 'str')
#dataset = pd.read_csv(f'{path}/TF predictivity.txt',sep='\t')
dataset = pd.read_csv('D:/TUT/Medical/biophysics/experiment/single cell sequencing/TF predictivity.txt',sep='\t')
seq1 = np.array(dataset.iloc[:,0].astype('str'),dtype = 'str')
#dataset2 = np.loadtxt('D:/TUT/Medical/biophysics/experiment/microarray sources.txt', delimiter ='\t',dtype = 'str')
#dataset2 = pd.read_csv(f'{path}/Partially reprogrammed cells Z-score.txt',sep='\t')
#seq2 = np.array(dataset2.iloc[:,0].astype('str'),dtype = 'str')
col = range(0,3)
row = range(1,np.size(dataset,1))
dataset3 = pd.read_csv('D:/TUT/Medical/biophysics/experiment/single cell sequencing/microarray sources.txt',sep='\t', usecols = col)
dataset1 = dataset.T.iloc[row,:]
cell_name = np.array(dataset1.index.astype('str'),dtype = 'str')
dataset4 = pd.read_csv(f'D:/TUT/Medical/biophysics/experiment/single cell sequencing/TF Z-score.txt',sep='\t')
seq3 = np.array(dataset4.iloc[:,0].astype('str'),dtype = 'str')
dataset4 = dataset4.T

class Walker:
    def __init__(self,*args,**kwargs):
        self.spin = kwargs['spin']
        self.nearest_neighbors = kwargs['nn']
        self.sys_dim = kwargs['dim']
        self.coords = kwargs['coords']

    def w_copy(self):
        return Walker(spin=self.spin.copy(),
                      nn=self.nearest_neighbors.copy(),
                      dim=self.sys_dim,
                      coords=self.coords.copy())
    

def Energy(Walkers):
    E = 0.0
    J = 4.0 # given in units of k_B
    for i in range(len(Walkers)):
        for k in range(len(Walkers[i].nearest_neighbors)):
            j = Walkers[i].nearest_neighbors[k]
            E += -J * Walkers[i].spin * Walkers[j].spin
    return E/2

def site_Energy(Walkers,s):
    E = 0.0
    J = 4.0 # given in units of k_B
#    i = len(Walkers)
    Walker = Walkers[s]
    for k in range(len(Walker.nearest_neighbors)):
        j = Walker.nearest_neighbors[k]                           
        E += -J*Walker.spin*Walkers[j].spin   
    return E

def ising(Nblocks,Niters,Walkers,beta):
    M = len(Walkers)
    Eb = np.zeros((Nblocks,))
    Accept=np.zeros((Nblocks,))
    AccCount=np.zeros((Nblocks,))
    N = np.zeros((Nblocks,))

    # Randomly change the orientation of a random cell in the lattice(upregulation/differentiation)
    # then record the energy changes and perform metropolis MC
    obs_interval = 1
    for i in range(Nblocks):
        EbCount = 1       
        NCount = 1
        for j in range(Niters):
            site = int(np.random.rand()*np.sqrt(M))
            
            s_old = 1.0*Walkers[site].spin
                
#            Walker[site] = 1.0*Walkers[site].spin
#            E_old = site_Energy(Walkers,Walkers[site])
            E_old = site_Energy(Walkers,site)
            # Choose randomly from the 2 spin choices:
            for l in range(site-1):
                Stemp = [2/(2**np.sqrt(M)-1)/(2**(site+1)-3), abs(1/(2**(site+1)-3)), (1+np.sqrt(M)-(site+1)-1+1/(2**(site+1)-3))/(2**np.sqrt(M)-1), 4/(2**np.sqrt(M)-1)*(1-1/(2**(site+1)-3)), 2**(site+1-l)/(2**np.sqrt(M)-1), (2**np.sqrt(M)-np.sqrt(M)-1-2**(site+1)+(site+1))/(2**np.sqrt(M)-1)]
            if j == np.sqrt(M):
                Stemp[2] = 1
            if j == 1:
                Stemp[2] = np.sqrt(M)/(2**(np.sqrt(M))-1)
            if  np.sort(Stemp)[-1] == 2**(site+1-l)/(2**np.sqrt(M)-1):
                N[i] = site+1
                s_new = np.sort(Stemp)[-1]
#            Walker[site] =np.sort(Stemp)[-1]*Walkers[site].spin 
                if Walker: del(Walker)
                E_new = s_new *site_Energy(Walkers, site)
            
                # Metropolis Monte Carlo:
                q_s = np.exp(-beta*(E_new-E_old))
                A_s = min(1.0, q_s)
                Ebdf = pd.DataFrame(Eb)
                Ebdf.fillna(0)
                Accdf = pd.DataFrame(Accept) 
                Accdf.fillna(0)
                Ndf = pd.DataFrame(N) 
                Ndf.fillna(0)
                Eb[i] = Ebdf.iloc[i]
                Accept[i] = Accdf.iloc[i]
                N[i] = Accdf.iloc[i]
                if (A_s > np.random.rand()):
                    Accept[i] += 1.0
                else:
                    Walkers[site].spin = 1.0*s_old
                
            
                AccCount[i] += 1.0
    
                if j % obs_interval == 0:
                    E_tot = Energy(Walkers)/np.sqrt(M) # energy per spin
                else:
                    E_tot = Eb[i-1] # energy per spin 
                Eb[i] += E_tot
                break
            
            s_new = np.sort(Stemp)[-1]
            
            if  j > obs_interval*10 and np.sort(Stemp)[-1] != 1 and np.sort(Stemp)[-1] == (2**np.sqrt(M)-np.sqrt(M)-1-2**(site+1)+(site+1))/(2**np.sqrt(M)-1):
                N[i] = 0
                E_tot = EbCount
                Accept[i] = 1
                Ebdf = pd.DataFrame(Eb)
                Ebdf.fillna(0)
                Accdf = pd.DataFrame(Accept) 
                Accdf.fillna(0)
                Ndf = pd.DataFrame(N) 
                Ndf.fillna(0)
                Eb[i] = Ebdf.iloc[i]
                Accept[i] = Accdf.iloc[i]
#                N[i] = Accdf.iloc[i]/ NCount
                Eb[i] += E_tot
                Accept[i] 
                AccCount[i] = 1.0
                break
            else:
                s_new = np.sort(Stemp)[-2]
            
#            Walker[site] =np.sort(Stemp)[-1]*Walkers[site].spin 
            E_new = s_new *site_Energy(Walkers,site)
            
                # Metropolis Monte Carlo:
            q_s = np.exp(-beta*(E_new-E_old))
            A_s = min(1.0, q_s)
            Ebdf = pd.DataFrame(Eb)
            Ebdf.fillna(0)
            Accdf = pd.DataFrame(Accept) 
            Accdf.fillna(0)
            Ndf = pd.DataFrame(N) 
            Ndf.fillna(0)
            Eb[i] = Ebdf.iloc[i]
            Accept[i] = Accdf.iloc[i]
            N[i] = Ndf.iloc[i]
            if (A_s > np.random.rand()):
                Accept[i] += 1.0
            else:
                Walkers[site].spin = 1.0*s_old
                N[i] += 1
                
            
            AccCount[i] += 1.0
    
            if j % obs_interval == 0:
                E_tot = Energy(Walkers)/np.sqrt(M) # energy per spin
            else:
                E_tot = Eb[i-1] # energy per spin 
            Eb[i] += E_tot
            EbCount += 1
            NCount += 1
            


        Ebdf = pd.DataFrame(Eb)
        Ebdf.fillna(0)
        Accdf = pd.DataFrame(Accept) 
        Accdf.fillna(0)
        Ndf = pd.DataFrame(N) 
        Ndf.fillna(0)
        Accept[i]
        N[i]
        Eb[i] = Ebdf.iloc[i]/ EbCount
        Accept[i] = Accdf.iloc[i]/ AccCount[i]
        N[i] = np.floor(Ndf.iloc[i]/ NCount)+1
        print('Block {0}/{1}'.format(i+1,Nblocks))
        print('    E   = {0:.5f}'.format(Eb[i]))
        print('    Acc = {0:.5f}'.format(Accept[i]))
        print('    Cellevel = {0:.5f}'.format(N[i]))

    return Walkers, Eb, Accept, N



def main():
    Walkers=[]

    dim = 2
    grid_side = np.size(dataset,1)
    grid_size = grid_side**dim
    
    # Ising model nearest neighbors only
    mapping = np.zeros((grid_side,grid_side),dtype=int) # mapping
    inv_map = [] # inverse mapping
    for i in range(grid_side):
        ii = 0
        for j in range(grid_side):
            mapping[i,j]=ii
            inv_map.append([i,j])
            ii += 1

    # Create a walker for each atom containing its coordinates, spin value and
    # its nearest neighbor mappings.
    for i in range(grid_side):
        for j in range(grid_side):
            j1=mapping[i,(j-1) % grid_side]
            j2=mapping[i,(j+1) % grid_side]
            i1=mapping[(i-1) % grid_side,j]
            i2=mapping[(i+1) % grid_side,j]
            Walkers.append(Walker(spin=1/6,
                                  nn=[j1,j2,i1,i2],
                                  dim=dim,
                                  coords = [i,j]))
 
    
    Nblocks = 200
    Niters = 1000
    eq = 600 # equilibration "time"
    T = 3.0
    beta = 1.0/T
    """
    Notice: Energy is measured in units of k_B, which is why
            beta = 1/T instead of 1/(k_B T)
    """
    # Perform classical monte carlo integration to get the ground state spin structure
    # and energy at given values.
    Walkers, Eb, Acc, N = ising(Nblocks,Niters,Walkers,beta)
    
    for c in range(1,int(max(N)+1)):
        plt.figure(),
        plt.subplot(221)
        plt.plot(abs(Eb[N==c])/max(abs(Eb[N==c])))
        plt.title([str(c),'cell level: '])
        plt.ylabel('expressions')
        plt.subplot(222)
        plt.plot(Accept[N==c])
        plt.ylabel('AcceptRatio')
        Etemp = Eb[N==c]
        ETtot = []
        ETmean = np.zeros(np.size(Etemp,0)-1)
        ETstd = np.zeros(np.size(Etemp,0)-1)
        for temp in range(1, np.size(Etemp,0)):
            ETtot.append(abs(Etemp[0:temp]))
            ETmean[temp-1] = np.mean(ETtot[temp-1])
            ETstd[temp-1] = 0.05*np.std(ETtot[temp-1])
        plt.subplot(223)
        plt.plot(ETstd/ETmean)
        plt.ylabel('Variance to energy ratio')
        plt.subplot(224)
        plt.plot(range(np.size(ETmean,0)),ETmean,'k')
        plt.fill_between(range(np.size(ETmean,0)),ETmean-ETstd[-1:], ETmean+ETstd[-1:])
        plt.ylabel('covergence')
        
        
#    plt.plot(Acc[Acc!=0])
#    plt.plot(N[N!=0])
#    plt.plot(Eb[Eb!=0])
#    plt.plot(Acc[Acc!=0])
#    plt.plot(N[N!=0])
#    Eb = Eb[eq:]
#    print('Ising total energy: {0:.5f} +/- {1:0.5f}'.format(np.mean(Eb), np.std(Eb)/np.sqrt(len(Eb))))
#    print('Variance to energy ratio: {0:.5f}'.format(abs(np.std(Eb)/np.mean(Eb)))) 
#    plt.show()

if __name__=="__main__":
    main()
