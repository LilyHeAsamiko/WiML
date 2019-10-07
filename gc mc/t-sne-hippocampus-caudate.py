# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:24:12 2019

@author: LilyHeAsamiko
"""
import numpy as np
import seaborn as sns; sns.set()

import sys; sys.path.append('D:\TUT\Medical\biophysics\experiment\genome\openTSNEPack')
from openTSNE.sklearn import TSNE
from openTSNE.callbacks import ErrorLogger

#from utils import utils
from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt

import pandas as pd
from sklearn import preprocessing
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
import h5py
import hdf5storage


#from numpy.random import Generate 
dataset = pd.read_csv(r'D:/TUT/Medical/biophysics/experiment/WiML/raw_data.tsv',sep='\t', header = 0)
dataset.to_hdf('D:/TUT/Medical/biophysics/experiment/WiML/raw_data.h5', 'table', append =  True)  
#dt = np.array(dataset, dtype = 'str')
#hdf5storage.write(dt, '.', 'raw_data.mat', matlab_compatible=True)

def cell2lb(cell):
    c = np.unique(cell)
    label = np.zeros((len(cell),1))
    for i in range(0, len(cell)):
        for j in range(0,len(c)):
            if cell[i] == c[j]:
                label[i] = j
    return label 

#dataset = pd.read_csv(r'D:/TUT/Medical/biophysics/experiment/WiML/raw_data.tsv',sep='\t', header = 0)


gene = dataset.iloc[:,1].astype('str')
cell = dataset.iloc[:,2].astype('str')
fpkm = dataset.iloc[:,3].astype('float32')
fpkm_norm = dataset.iloc[:,4].astype('float32')
Label = cell2lb(cell)
x_train = np.zeros((len(cell),2))
x_train[:,0] = fpkm_norm
x_train[:,1] = np.array(cell2lb(gene)).ravel()

lb = preprocessing.LabelBinarizer()
lb.fit(Label)
lb_b = lb.transform(Label)
y_train = lb_b

mu = np.mean(x_train)
sigma = np.std(x_train)
dataset = np.zeros((75, 25))
for i in range(np.shape(dataset)[1]):
    dataset[:,i] = np.random.normal(mu, sigma, 75)

#count0, bin0, ignored = plt.hist(x_train, col) 
#x_train, x_test, y_train, y_test = train_test_split(x_train, y_train, test_size=33/75, random_state=42)

#(x_train, y_train), (x_test, y_test) = mnist.load_data()
#x_train = x_train.reshape(60000, 784).astype('float64') / 255
#x_test  =  x_test.reshape(10000, 784).astype('float64') / 255


tsne = []

tsne = TSNE(
#    perplexity=30, 5-50
    perplexity = 10,
    metric="euclidean",
    callbacks=ErrorLogger(),
    n_jobs=20, #minimum 250(minimum iterations)
    random_state=42, #10-1000 might increase to avoid cost function being stuck in a bad local minimum
#    init = 'pca'
)

#col = np.repeat('0000000',np.size(lb_b,1))
#for i in range(np.size(lb_b,1)):
#    R = i
#    G = i+20
#    B = i+40
#    col[i] ='#{:02X}{:02X}{:02X}'.format(R, G, B)
#color
#col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
#                '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']),45)

col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#c6cee3','#1b78b4','#c2df8a','#34d02c']),3)

embedding_train = tsne.fit(x_train)
fig = plt.figure()
plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1], c=col, s=1)
plt.xlim(-300, 200)
plt.ylim(min(embedding_train.embedding_[:,1]), 300)
#plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1]), # c=col[Label.astype('uint32')], s=1)
plt.tight_layout()
plt.savefig('training tsne.pdf', dpi=200)
           
dataset2 = pd.read_csv(r'D:/TUT/Medical/biophysics/experiment/WiML/MAdata_DGinMCout.txt',sep='\t', header = 0)
cell = np.repeat([1,2,3],25).ravel();
fpkm_post = np.array(dataset2.astype('float32')).ravel()
x_train = np.zeros((len(cell),2))
x_train[:,0] = fpkm_post
x_train[:,1] = np.array(cell2lb(gene)).ravel()
Label = cell

lb = preprocessing.LabelBinarizer()
lb.fit(Label)
lb_b = lb.transform(Label)
y_train = lb_b

x_train, x_test, y_train, y_test = train_test_split(x_train, y_train, test_size=.25, random_state=42)

#(x_train, y_train), (x_test, y_test) = mnist.load_data()
#x_train = x_train.reshape(60000, 784).astype('float64') / 255
#x_test  =  x_test.reshape(10000, 784).astype('float64') / 255

tsne = []

tsne = TSNE(
#    perplexity=30, 5-50
    perplexity = 10,
    metric="euclidean",
    callbacks=ErrorLogger(),
    n_jobs=20, #minimum 250(minimum iterations)
    random_state=42, #10-1000 might increase to avoid cost function being stuck in a bad local minimum
#    init = 'pca'
)

#col = np.repeat('0000000',np.size(lb_b,1))
#for i in range(np.size(lb_b,1)):
#    R = i
#    G = i+20
#    B = i+40
#    col[i] ='#{:02X}{:02X}{:02X}'.format(R, G, B)
#color
#col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
#                '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']),45)

col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#c6cee3','#1b78b4','#c2df8a','#34d02c']),4)

embedding_train = tsne.fit(x_train)

fig = plt.figure()
plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1], c=col, s=1)
plt.xlim(-300, 200)
plt.ylim(min(embedding_train.embedding_[:,1]), 300)
#plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1]), # c=col[Label.astype('uint32')], s=1)
plt.tight_layout()

X = np.concatenate((x_train, x_test))
y = np.concatenate((y_train, y_test))

# Do PCA and keep 50 dimensions
X = X - X.mean(axis=0)
U, s, V = np.linalg.svd(X, full_matrices=False)
X50 = np.dot(U, np.diag(s))[:,:50]

#initialization
PCAinit = X50[:,:2] / np.std(X50[:,0]) 

col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#c6cee3','#1b78b4','#c2df8a','#34d02c','#17d01a']),5)

embedding_test = tsne.fit(PCAinit)
fig = plt.figure()
plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1], c=col, s=1)
plt.xlim(-300, 200)
plt.ylim(min(embedding_train.embedding_[:,1]), 300)
#plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1]), # c=col[Label.astype('uint32')], s=1)
plt.tight_layout()

plt.savefig('test tsne after PCA.pdf', dpi=200)

from sklearn.cluster import KMeans
Col = np.shape(dataset)[1]
Col = 4
#    def()
for col in range(2, Col):

    
#    iters = 10
    iters = 2
    KL = np.zeros((len(dataset),len(dataset), col, iters))
    SS = np.zeros((len(dataset),len(dataset), col, iters)) 
    std_KL = np.zeros((len(dataset),len(dataset), col, iters))
    std_SS = np.zeros((len(dataset),len(dataset), col, iters))  
    Perp = np.zeros((len(dataset),len(dataset), col, iters))
        #    delta_P = 1000*np.zeros((len(dataset), col))
        #    delta_P = np.repeat(1000, len(dataset), axis = 0)
    Cluster_Perp = np.zeros((len(dataset), col, iters))
    Cluster = np.zeros((np.shape(dataset)[0], np.shape(dataset)[1], col, iters))
    Np = np.zeros((col, iters))
    #Npi = np.zeros((np.shape(dataset)[1], 2, 1, 1))
    P = np.zeros((np.shape(dataset)[0], np.shape(dataset)[0], np.shape(dataset)[1], 2, iters))
    #Trans = np.zeros((np.shape(dataset)[0], np.shape(dataset)[0], np.shape(dataset)[1], 2, 1, 1))
    Bin = np.zeros((np.shape(dataset)[0], np.shape(dataset)[0], col, iters)) 
    delta_P = np.zeros((np.shape(dataset)[0], col, iters))
#    if col == 1:
#        mu = np.mean(x_train)
#        sigma = np.std(x_train)       
#    else:
    for step in range(iters):
#        print(['step', step])
        count0, bin0, ignored = plt.hist(dataset.T, col, align = 'mid') 
    #    dt1_df['cluster_KNN'] = cluster
    #    Cluster[col-1,] = cluster
    #    ii = 0
    #    for i in range(len(x_train)):
    #    bin0 = bin0[1:col+1,]
        bin0 = (bin0[0:col,] + bin0[1:col+1,])/2
        d = np.array(count0)*np.array(bin0)/np.shape(dataset)[1]
        S = 0
        E = 0.00001
    #    sigma = 0
    
        for s in range(np.shape(d)[0]):
            for t in range(np.shape(d)[0]):
                if s != t:
                    S += 1/(1+(d[s]-d[t])**2)                                
    #                print(S)
        for s in range(np.shape(dataset)[0]):
            for t in range(np.shape(dataset)[0]):
                if s != t:
                    sigma = np.std(dataset[t,])
                    E += np.exp(-(dataset[s]-dataset[t])**2/(2*sigma**2+0.00001))
    #                print(E)
        for i in range(len(dataset)):
            print(['i',i])
            for j in range(len(dataset)):
                print(['step:',step])
                print(['j',j])
                print(['i',i])
                p_i = np.zeros((1,np.shape(dataset)[1]))
                p_j = np.zeros((1,np.shape(dataset)[1]))
                perp = np.zeros((1,col)) 
                if j != i:
                    sigma = np.std(dataset[i,])   
                    sigma1 = np.std(dataset[j,]) 
                    kl =  np.zeros((1,col))
                    H =  np.zeros((1,col)) 
                    q = 1/(1+(d[i,]-d[j,])**2)/S
#                    p_i = np.exp(-(dataset[i,]-dataset[j,])**2/(2*sigma1**2+0.0001))/E
#                    p_j = np.exp(-(dataset[j,]-dataset[i,])**2/(2*sigma**2+0.0001))/E
                    p_i = np.exp(-(dataset[i,]-dataset[j,])**2/(2*sigma1**2+0.0001))
                    p_j = np.exp(-(dataset[j,]-dataset[i,])**2/(2*sigma**2+0.0001))
                    print(p_i)
#                    print(p_j)
                    p = (p_i + p_j)/(2*np.shape(dataset)[1])
                    P[i, j, :, 0,  step] = np.array(p_i)
                    P[i, j, :, 1,  step] = np.array(p_j)
                    countp, binp, ignored = plt.hist(p, col)
                    countE, binE, ignored = plt.hist(p*np.log10(E), col)
                    countEE, binEE, ignored = plt.hist(E, col)
    #                        p = (p+np.exp(-(x_train[j]-x_train[i])**2/(2*sigma1**2+0.0001))/E)/(2*len(x_train))
                    binp = binp[1:col+1, ] 
                    binE = binE[1:col+1, ]  
                    binEE = binEE[1:col+1, ] 
#                    print(binp)
#                    print(binE)
#                    print(binEE)
#                    Bin[i, j, :, step] = binp
#                    Trans[i, j, :, col-2, step] = p/binp
#                    kl = binp*np.log10(binp/q) - binp*np.log10(E)
                    kl = binp*np.log10(binp/q) - binE
                    H = -sigma1*binp*np.log2(binp) + sigma1*binp*np.log2(binEE) 
                    perp = 2**H # somehow has to be defined first
    #                print(perp)
#                    Perp[i, j, :, col-2, step] = perp
                    SS[i, j, :, step] = 2*kl + np.log10(len(dataset))*perp/(len(dataset))
                    KL[i, j, :, step] = kl
                    Perp[i, j, :, step] = np.nan_to_num(perp)
#                    Perp[0, 0, step, ] = Perp[0, 1, ]
                    SS = np.nan_to_num(SS)
                    KL = np.nan_to_num(KL)
    #                std_SS[i, j]= 3*np.std(SS[np.shape(SS)!=0])
    #                std_KL[i, j]= 3*np.std(KL[np.shape(KL)!=0])                
                    std_SS[i, j, :, step]= 3*np.std(SS[i, 1:j, :,  step])
                    std_KL[i, j, :, step]= 3*np.std(KL[i, 1:j, :, step])
                    delta_P[: , :,  step] = abs(Perp[0, :, :, step]-Perp[:, 0, :, step])
#                    print(['delta_P', delta_P])
    
                    if j <= i:
                        for c in range(1, col):
                            print(['c',c])
                            delta_Perp = abs(Perp[j, :, c-1, step] - Perp[:, j, c-1, step])
#                            print(['dPerp=', sum(delta_Perp)])
#                            print(['dP=', sum(delta_P[:, c-1, step])])
#                            print( sum(delta_Perp) < sum(delta_P[:, c-1, step]) and i*j!=0)
                            if sum(delta_Perp) < sum(delta_P[:, c-1, step]) and i*j!=0:
                                print(['j =',j])
                                temp_equation = Np[c-1, step] == j
                                if np.size(temp_equation) == 1: 
                                    if Np[c-1, step] == j:
                                        break                                
                                else:
                                    if sum(Np[c-1, step] == j) == 0: 
         #                               print(j)
         #                               print(Np[c-1])
                                        delta_P[:, c-1, step] = delta_Perp
                                        Cluster_Perp[:, c-1, step] = (Perp[j, :, c-1, step] + Perp[:, j, c-1, step])/2  
#                                        print(j)
                                        Np[c-1, step] = j
                                        Sigma1 = np.log(Cluster_Perp[j, c-1, step]+0.000001, np.array(-Bin[j, j, c-1, step]))/(Bin[j, j, c-1, step]+0.00001)
        #                                P1 = Trans[j, j, :, c-1]* Bin[j, j, c-1]
                                        if np.size(Sigma1) == 1: 
                                            Data = -np.log(np.e, P[j, j, :, 0, step].T * E)*2*Sigma1**2 + dataset[j, ]
                                        else:
                                            Data = -np.log(np.e, P[j, j, :, 0, step].T * E)*2*Sigma1[c-1]**2 + dataset[j, ]
        #                                delta_NPerp = dataset[j, :]-Cluster_Perp[:, c-1]
        #                                Ind = np.argmin(delta_NPerp)
#                                        Cluster[j, :, c-1, step] = Data
#                                        print(Cluster[j, :, c-1, step])
    #                                    Npi[:, c-1, col-2, step] = np.argmin(dataset[j,] - Data[0])

        Step = np.zeros((np.shape(Np)[0], np.shape(Np)[1]))
        Step = np.array(Step, dtype = 'int32') 
        temp = SS[0, :, 0, Step[0,0]]
        col_S = 0
        i_S = 0
        for c in range(col):
            for i in np.array(Np[c,], dtype = int):
                i = np.int32(i)
                Step[c,i] = np.int32(np.argmin(KL[i,i,c,:]))
                Step = np.array(Step, dtype = 'int32')                
                plt.figure(),
                plt.subplot(211)
                plt.plot(np.linspace(0, np.size(KL[i, :, c, Step[c,i]]), np.size(KL[i, :, c, Step[c,i]])), np.array(KL[i, :, c, Step[c,i]]), 'r')
                plt.fill_between(np.linspace(0, np.size(KL[i, :, c, Step[c,i]]), np.size(KL[i, :, c, Step[c,i]])), np.array(KL[i, :, c, Step[c,i]] - std_KL[i, :, c, Step[c,i]]), np.array(KL[i, :, c, Step[c,i]] + std_KL[i, :, c, Step[c,i]]), color = 'blue', alpha = 0.05)
                plt.axvline(x = i, ymin =  min(KL[i, :, c, Step[c,i]] - std_KL[i, :, c, Step[c,i]])-0.02, ymax = max(KL[i, :, c, Step[c,i]]+ std_KL[i, :, c, Step[c,i]])+0.02)
                plt.title('KL')
                plt.subplot(212)
                plt.plot(np.linspace(0, np.size(SS[i, :, c, Step[c,i]]), np.size(KL[i, :, c, Step[c,i]])), np.array(KL[i, :, c, Step[c,i]]), 'r')
                plt.fill_between(np.linspace(0, np.size(SS[i, :, c, Step[c,i]]), np.size(SS[i ,:, c, Step[c,i]])), np.array(SS[i, :, c, Step[c,i]] - std_SS[i, :, c, Step[c,i]]), np.array(SS[i, :, c, Step[c,i]] + std_SS[i, :, c, Step[c,i]]), color = 'blue', alpha = 0.05)
                plt.axvline(x = i, ymin =  min(SS[i, :, c, Step[c,i]]- std_SS[i, :, c, Step[c,i]])-0.02, ymax = max(SS[i, :, c, Step[c,i]] + std_SS[i, :, c, Step[c,i]])+0.02)
                plt.title('S(perplex)')
                plt.savefig('KL_S(perplex)_%d.png' % np.int32(col*(c-1)+i))
                if sum(temp) > sum(SS[i, :, c, Step[c,i]]):
                    temp = SS[i, :, c, Step[c,i]]
                    col_S = c
                    i_S = i
#                interv = np.int32((np.shape(dataset)[1]+1)/col)
                temp_col = dataset[:,abs(np.array(dataset[i,]) - bin0[c])< np.mean(abs(np.array(dataset[i,]) - bin0[c]))]
   #                Cluster[i, -1-(c+1)*interv:-1-c*interv, c, Step[c,i]] = temp_col[0:len(Cluster[i, -1-(c+1)*interv:-1-c*interv, c, Step[c,i]])]                
                Cluster[i, abs(np.array(dataset[i,]) - bin0[c])< np.mean(abs(np.array(dataset[i,]) - bin0[c])), c, Step[c,i]] = abs(np.array(dataset[i,]) - bin0[c])[ abs(np.array(dataset[i,]) - bin0[c])< np.mean(abs(np.array(dataset[i,]) - bin0[c]))]
    #                plt.savefig('D:/TUT/Medical/biophysics/experiment/WiML/KL and S(perplex) with %d th cluster and center at %d th RNA.png ' % (c, i), dpi = 200)
#                plt.savefig('D:/TUT/Medical/biophysics/experiment/WiML/KL and S(perplex)_%d.png ' % np.int32(col*c+i), dpi = 200
#        C = abs(Cluster[i_S,i_S,:,Step[col_S,i_S]])/np.max(abs(Cluster[i_S,i_S,:,Step[col_S,i_S]]))
        C = abs(Cluster[i_S, :, :, Step[col_S,i_S]]) 
        C[C==0] = 100
        C = np.argmin(C, axis =1)
        C = C/max(C)
# [XX,YY] = np.sort(np.meshgrid(np.array(dataset[i_S,], dtype = int),np.array(dataset[:,i_S], dtype = int).T))
        plt.figure(),
        plt.subplot(111)
        plt.scatter(dataset[i_S,], dataset[i_S+1,], c = C, cmap = 'coolwarm')
        plt.show()
        plt.savefig('custer_%d.png' % np.int32(col*(col_S-1)+i_S))


            
    Np = np.array(Np, dtype = 'int32')
    KL[KL == 0]=0.5
    SS[SS == 0]=0.5
    KL_col = np.argmin(KL, axis = 2)
    SS_col = np.argmin(SS, axis = 2)
    KL_col[KL_col == 0]=np.max(KL_col)
    SS_col[SS_col == 0]=np.max(SS_col)
    plt.figure,
    plt.subplot(211)
    plt.plot(np.linspace(KL_col, np.size(KL[Np[i_S, np.min(KL_col)], 0:16, np.min(Step[KL_col,])], np.size(KL[Np[i_S, np.min(KL_col)], 0:16, np.min(Step[KL_col,])]), KL[Np[i_S, np.min(KL_col)], 0:16, np.min(Step[KL_col,])])    
    plt.subplot(212)
    plt.plot(np.linspace(SS_col, np.size(SS[Np[i_S, np.min(SS_col)], 0:16, np.min(Step[SS_col,])], np.size(KL[Np[i_S, np.min(KL_col)], 0:16, np.min(Step[KL_col,])]), SS[Np[i_S, np.min(SS_col)], 0:16, np.min(Step[SS_col,])])
    
    
#    C = abs(Cluster[:,:,0,0])/np.max(abs(Cluster[:,:,0,0]))
#    [XX,YY] = np.sort(np.meshgrid(np.array(dataset, dtype = int),np.array(dataset, dtype = int).T))
#    fig = plt.figure()
#    plt.imshow(np.meshgrid(np.array(C), np.array(C).T)[0])
#    plt.show()
#    plt.savefig('cell evolution_colormap.pdf',dpi=200)
#    
#    plt.figure()    
#    plt.subplot(111, projection='3d')
##    plt.scatter(XX,YY, np.meshgrid(np.array(C), np.array(C).T)[0])
#    plt.scatter(XX, YY, C)
#    plt.title('cell_evolution(cell1) scatter on x_z_y')
##    plt.show()
#    plt.figure
#    plt.subplot(111, projection='3d')
#    x,y = np.meshgrid(dataset[:,0].astype('float16'), np.linspace(min(dataset[0,:].astype('float16')),max(dataset[0,:].astype('float16')),1000))
#    z = Cluster[:,:,0]
#    for c, zlow, zhigh in [(1, min(z[0,:]), max(z[0,:])),(2, min(z[1,:]),max(z[1,:])),(3, min(z[2,:]),max(z[2,:])),(4, min(z[3,:]),max(z[3,:]))]:
#        plt.scatter(x[0,:].astype('float16'), y[0, :],  c, cmap='prism')  # plot points with cluster dependent colors
#    plt.title('cell_evolution(4RNA) scatter on x_z_y')
#    plt.show()
#    fig.savefig('cell evolution_3D(4RNA).pdf',dpi=200)    
    

dataset = pd.read_csv(r'D:/TUT/Medical/biophysics/experiment/WiML/raw_data_DGinMCout.tsv',sep='\t', header = 0)
gene = dataset.iloc[:,1].astype('str')
cell = dataset.iloc[:,2].astype('str')
fpkm = dataset.iloc[:,3].astype('float32')
fpkm_norm = dataset.iloc[:,4].astype('float32')
Label = cell2lb(cell)
x_train = np.zeros((len(cell),2))
x_train[:,0] = fpkm_norm
x_train[:,1] = np.array(cell2lb(gene)).ravel()

lb = preprocessing.LabelBinarizer()
lb.fit(Label)
lb_b = lb.transform(Label)
y_train = lb_b

x_train, x_test, y_train, y_test = train_test_split(x_train, y_train, test_size=33/75, random_state=42)

#(x_train, y_train), (x_test, y_test) = mnist.load_data()
#x_train = x_train.reshape(60000, 784).astype('float64') / 255
#x_test  =  x_test.reshape(10000, 784).astype('float64') / 255

tsne = []

tsne = TSNE(
#    perplexity=30, 5-50
    perplexity = 10,
    metric="euclidean",
    callbacks=ErrorLogger(),
    n_jobs=20, #minimum 250(minimum iterations)
    random_state=42, #10-1000 might increase to avoid cost function being stuck in a bad local minimum
#    init = 'pca'
)

#col = np.repeat('0000000',np.size(lb_b,1))
#for i in range(np.size(lb_b,1)):
#    R = i
#    G = i+20
#    B = i+40
#    col[i] ='#{:02X}{:02X}{:02X}'.format(R, G, B)
#color
#col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
#                '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']),45)

col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a']), 14)#,'#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#c6cee3','#1b78b4','#c2df8a','#34d02c']),3)

embedding_train = tsne.fit(x_train)
fig = plt.figure()
plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1], c=col, s=1)
plt.xlim(-300, 200)
plt.ylim(min(embedding_train.embedding_[:,1]), 300)
#plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1]), # c=col[Label.astype('uint32')], s=1)
plt.tight_layout()
plt.savefig('training tsne.pdf', dpi=200)
           
dataset2 = pd.read_csv(r'D:/TUT/Medical/biophysics/experiment/WiML/MAdata_DGinMCout.txt',sep='\t', header = 0)
cell = np.repeat([1,2,3],25).ravel();
fpkm_post = np.array(dataset2.astype('float32')).ravel()
x_train = np.zeros((len(cell),2))
x_train[:,0] = fpkm_post
x_train[:,1] = np.array(cell2lb(gene)).ravel()
Label = cell

lb = preprocessing.LabelBinarizer()
lb.fit(Label)
lb_b = lb.transform(Label)
y_train = lb_b

x_train, x_test, y_train, y_test = train_test_split(x_train, y_train, test_size=33/75, random_state=42)

#(x_train, y_train), (x_test, y_test) = mnist.load_data()
#x_train = x_train.reshape(60000, 784).astype('float64') / 255
#x_test  =  x_test.reshape(10000, 784).astype('float64') / 255

tsne = []

tsne = TSNE(
#    perplexity=30, 5-50
    perplexity = 10,
    metric="euclidean",
    callbacks=ErrorLogger(),
    n_jobs=20, #minimum 250(minimum iterations)
    random_state=42, #10-1000 might increase to avoid cost function being stuck in a bad local minimum
#    init = 'pca'
)

#col = np.repeat('0000000',np.size(lb_b,1))
#for i in range(np.size(lb_b,1)):
#    R = i
#    G = i+20
#    B = i+40
#    col[i] ='#{:02X}{:02X}{:02X}'.format(R, G, B)
#color
#col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
#                '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']),45)

col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a']),14)#,'#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#c6cee3','#1b78b4','#c2df8a','#34d02c']),4)

embedding_train = tsne.fit(x_train)

fig = plt.figure()
plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1], c=col, s=1)
plt.xlim(-300, 200)
plt.ylim(min(embedding_train.embedding_[:,1]), 300)
#plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1]), # c=col[Label.astype('uint32')], s=1)
plt.tight_layout()









 
X = np.concatenate((x_train, x_test))
y = np.concatenate((y_train, y_test))

# Do PCA and keep 50 dimensions
X = X - X.mean(axis=0)
U, s, V = np.linalg.svd(X, full_matrices=False)
X50 = np.dot(U, np.diag(s))[:,:50]

#initialization
PCAinit = X50[:,:2] / np.std(X50[:,0]) 

col = np.repeat(np.array(['#a6cee3','#1f78b4','#b2df8a']),14)#,'#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#c6cee3','#1b78b4','#c2df8a','#34d02c','#17d01a']),5)

embedding_test = tsne.fit(PCAinit)
fig = plt.figure()
plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1], c=col, s=1)
plt.xlim(-300, 200)
plt.ylim(min(embedding_train.embedding_[:,1]), 300)
#plt.scatter(embedding_train.embedding_[:,0], embedding_train.embedding_[:,1]), # c=col[Label.astype('uint32')], s=1)
plt.tight_layout()

plt.savefig('test tsne after PCA.pdf', dpi=200)

