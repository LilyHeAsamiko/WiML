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

def prune(P0, j, L, Seq, k, LDNA, l, c, fb):
    Pf_ = np.array(list([0]*len(x)))
    Pb_ = np.array(list([0]*len(x)))
    L = L-1;
    if fb == 'f' :
        if k - 1 ==  k - L:
            nfn_ = Seq[k-1].count(LDNA[j])
        else:
            nfn_ = Seq[k-1:k-L].count(LDNA[j])
        if nfn_ == 0:
            nfn_ = Seq[k-5:k+5-1].count(LDNA[j])/len(LDNA[j])**(l+2)+0.00001 
        if k ==  k - L:
            nfd_ = Seq[k].count(LDNA[j])
        else:
            nfd_ = Seq[k-L:k].count(LDNA(j))
        if nfd_ == 0:
            nfd_ = Seq[k-5:k+5-1].count(LDNA[j])/len(LDNA[j])**(l+2)+0.00001 
        if nfd_ < 0.001:
            if nfn_ == 0:
                Pf_ = 0
            else:
                Pf_ = 1
        else:
            Pf_ = nfn_/nfd_   
        
        KL = sum(Pf_*np.log(2, P0/Pf))   
        if KL <= c:
            Pf_= 0
        P = Pf_
            
        
    else:
        if k - 1 ==  k - Lb:
            nbn_ = Seq[k-1].count(dDNA[j])
        else:
            nbn_ = Seq[k-L:k-1].count(dDNA(j))
        if nbn_ == 0:
            nbn_ = Seq[k-5:k+5-1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
        if k - 1 ==  k - L:
            nbn_ = Seq[k-1].count(dDNA[j])
        else:
            nbd_ = Seq[k-L:k].count(dDNA(j))
        if nbd_ == 0:
            nbd_ = Seq[k-5:k+5-1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
        if nfd_ < 0.001:
            if nbn_ == 0:
                Pb_ = 0
            else:
                Pb_ = 1
        else:
            Pb_ = nbn_/nbd_   
        
        KL = sum(Pb_*np.log(2, P0/Pb))  
        psi = len(LDNA)**(l+1)/nbn_
        while L > 1 & KL <= c*psi:
            L = L - 1
            Pb_ = prune(Pb_, L, Seq, k, dDNA, l, c, 'b')
        P = Pb_
    return P
#path = r'D:/Introduction to R/HMM//'
file = r'P:/GCMC/GCMC/STXBP1pheno.txt'
#file = r'D:/Introduction to R/HMM/STXBP1pheno.txt'
dataset = pd.read_csv(file, sep = ',')
Seq = str(dataset.x[1])
seq = list(dataset.x[1])
x = seq[0:99] #example 100
        

dDNA = ('A','T','C','G')
Lf = 1
Lb = 5# binding site
c = 0.65
t = list([]*len(dDNA))
parents = list([''])
childsN =list([])

Pf = np.array(list([0.0]*len(x)))
Pb = np.array(list([0.0]*len(x)))
nfn = np.array(list([0.0]*len(dDNA)))
nfd = np.array(list([0.0]*len(dDNA)))
nbn = np.array(list([0.0]*len(dDNA)))
nbd = np.array(list([0.0]*len(dDNA)))
PF = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))
PB = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))

pred_1 =np.array(list(['']*len(x)))
accA = np.array(list([0.0]*len(x)))
accC = np.array(list([0.0]*len(x)))
accT = np.array(list([0.0]*len(x)))
accG = np.array(list([0.0]*len(x)))
acc = np.array([accA,accC,accT,accG])
tfbsA = np.array(list([0.0]*len(x)))
tfbsT = np.array(list([0.0]*len(x)))
tfbsC = np.array(list([0.0]*len(x)))
tfbsG = np.array(list([0.0]*len(x))) 

#find tfbs with 10 box(fragment with significantly higher accuracy)
for i in np.arange(1, len(x)):  
    for j in range(len(dDNA)):
        if i - 1 ==  i - Lf:
            nfn[j] = Seq[i-1].count(dDNA[j])
        else:
            #[]: upper bound excluded
            nfn[j] = Seq[i-Lf:i-1+1].count(dDNA(j))
        if nfn[j] == 0:
            nfn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(i+1)+0.00001 
        if i ==  i - Lf:
            nfd[j] = Seq[i].count(dDNA[j]) 
        else:
            nfd[j] = Seq[i-Lf:i+1].count(dDNA[j])
        if nfd[j] == 0:
            nfd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(i+1)+0.00001 
        if nfd[j] < 0.001:
            if nfn[j] == 0:
                Pf[j] = 0
            else:
                Pf[j] = 1
        else:
            Pf[j] = nfn[j]/nfd[j]        
        if  x[i] == dDNA[0]:
            PF[0,j] = Pf[j]
    #            PB[0,j] = max(Pb[j], Pb_[j])               
        if  x[i] == dDNA[1]: 
            PF[1,j] = Pf[j]
    #            PB[1,j] = max(Pb[j], Pb_[j])
        if  x[i] == dDNA[2]: 
            PF[2,j] = Pf[j]
    #            PB[2,j] = max(Pb[j], Pb_[j])
        if  x[i] == dDNA[3]: 
            PF[3,j] = Pf[j]
    #            PB[3,j] = max(Pb[j], Pb_[j])     
        if j == 3:
            nfn = np.array(list([0]*len(dDNA)))
            nfd = np.array(list([0]*len(dDNA)))
#            nbn_ = np.array(list([0]*len(LDNA)))
#            nbd_ = np.array(list([0]*len(LDNA)))
#            nbn = np.array(list([0]*len(LDNA)))
#            nbd = np.array(list([0]*len(LDNA))) 
    P0 = np.dot([1,1,1,1], PF)
    pred_1[i] = dDNA[np.argmax(P0)]
    accA[i] = sum((pred_1[0:i+1]==x[0:i+1]) * (pred_1[0:i+1] == np.repeat(dDNA[0],i+1)))/(i+1)
    accT[i] = sum((pred_1[0:i+1]==x[0:i+1]) * (pred_1[0:i+1] == np.repeat(dDNA[1],i+1)))/(i+1)
    accC[i] = sum((pred_1[0:i+1]==x[0:i+1]) * (pred_1[0:i+1] == np.repeat(dDNA[2],i+1)))/(i+1)
    accG[i] = sum((pred_1[0:i+1]==x[0:i+1]) * (pred_1[0:i+1] == np.repeat(dDNA[3],i+1)))/(i+1)        
acc = np.array([accA,accC,accT,accG])
tfbsA[(accA>np.mean(accA))*(accA>=accC)*(accA>=accT)*(accA>=accG)*(pred_1==dDNA[0])] = 1
tfbsT[(accT>np.mean(accT))*(accT>=accA)*(accT>=accC)*(accT>=accG)*(pred_1==dDNA[1])] = 1
tfbsC[(accC>np.mean(accC))*(accC>=accA)*(accC>=accT)*(accC>=accG)*(pred_1==dDNA[2])] = 1
tfbsG[(accG>np.mean(accG))*(accG>=accA)*(accG>=accC)*(accG>=accT)*(pred_1==dDNA[3])] = 1
tfbs = np.array([tfbsA,tfbsC,tfbsT,tfbsG]) 

Pf = np.array(list([0.0]*len(x)))
Pb = np.array(list([0.0]*len(x)))
nfn = np.array(list([0.0]*len(dDNA)))
nfd = np.array(list([0.0]*len(dDNA)))
nbn = np.array(list([0.0]*len(dDNA)))
nbd = np.array(list([0.0]*len(dDNA)))
PF = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))
PB = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))

#fix L, 10 datapoints, homogenious VOM(5,0.65) tree
for i in np.arange(1, len(x)):
    l = 0
    if(sum(tfbs[:,i] == 0)>0 & i >= Lb):
        for j in range(len(dDNA)):
            if i - 1 ==  i - Lb:
                nbn[j] = Seq[i-1].count(dDNA[j])
            else:
                nbn[j] = Seq[i-Lb:i-1+1].count(dDNA(j))
            if nbn[j] == 0:
                nbn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
            if i ==  i - Lb:
                nbd[j] = Seq[i].count(dDNA[j])
            else:       
                nbd[j] = Seq[i-Lb:i+1].count(dDNA(j))
            if nbd[j] == 0:
                nbd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
            if nbd[j] < 0.001:
                if nbn[j] == 0:
                    Pb[j] = 0
                else:
                    Pb[j] = 1
            else:
                Pb[j] = nbn[j]/nbd[j]            
            if  x[i] == dDNA[0]:
                PB[0,j] = Pb[j]
            #            PB[0,j] = max(Pb[j], Pb_[j])               
            if  x[i] == dDNA[1]: 
                PB[1,j] = Pb[j]
            #            PB[1,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[2]: 
                PB[2,j] = Pb[j]
            #            PB[2,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[3]: 
                PB[3,j] = Pb[j] 
            #            PB[3,j] = max(Pb[j], Pb_[j])  
            if j == 3:
                nbn = np.array(list([0]*len(dDNA)))
                nbd = np.array(list([0]*len(dDNA)))
    #            nbn_ = np.array(list([0]*len(LDNA)))
    #            nbd_ = np.array(list([0]*len(LDNA)))
    #            nbn = np.array(list([0]*len(LDNA)))
    #            nbd = np.array(list([0]*len(LDNA)))             
    #    t.append(np.dot(P, [1,1,1,1]))
        P0 = np.dot([1,1,1,1], PB)
        for j in range(len(dDNA)):    
            if i - 1 ==  i - Lb +1:
                nbn[j] = Seq[i-1].count(dDNA[j])
            else:
                nbn[j] = Seq[i-Lb+1:i-1+1].count(dDNA(j))
            if nbn[j] == 0:
                nbn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
            if i ==  i - Lb:
                nbd[j] = Seq[i].count(dDNA[j])
            else:       
                nbd[j] = Seq[i-Lb+1:i+1].count(dDNA(j))
            if nbd[j] == 0:
                nbd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
            if nbd[j] < 0.001:
                if nbn[j] == 0:
                    Pb[j] = 0
                else:
                    Pb[j] = 1
            else:
                Pb[j] = nbn[j]/nbd[j]           
            psi = len(dDNA)**(l+1)/nbn[j]
            KL = sum(Pb*np.log(2, P0/Pb))
            if KL <= c*psi:
                Pb[j] = prune(Pb[j], j, Lb, Seq, i, dDNA, l, c, 'b')        
            if x[i] == dDNA[0]:
                PB[0,j] = Pb[j]
            #            PB[0,j] = max(Pb[j], Pb_[j])               
            if  x[i] == dDNA[1]: 
                PB[1,j] = Pb[j]
            #            PB[1,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[2]: 
                PB[2,j] = Pb[j]
            #            PB[2,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[3]: 
                PB[3,j] = Pb[j]
            #            PB[3,j] = max(Pb[j], Pb_[j])     
            if j == 3:
                nfn = np.array(list([0]*len(dDNA)))
                nfd = np.array(list([0]*len(dDNA)))

            
    else:
        for j in range(len(dDNA)):
            if i - 1 ==  i - Lf:
                nfn[j] = Seq[i-1].count(dDNA[j])
            else:
                nfn[j] = Seq[i-Lf:i-1+1].count(dDNA(j))
            if nfn[j] == 0:
                nfn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
            if i  ==  i - Lf:
                nfd[j] = Seq[i].count(dDNA[j])
            else:
                nfd[j] = Seq[i-Lf:i+1].count(dDNA(j))
            if nfd[j] == 0:
                nfd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
            if nfd[j] < 0.001:
                if nfn[j] == 0:
                    Pf[j] = 0
                else:
                    Pf[j] = 1
            else:
                Pf[j] = nfn[j]/nfd[j]            
            if  x[i] == dDNA[0]:
                PF[0,j] = Pf[j]
    #            PB[0,j] = max(Pb[j], Pb_[j])               
            if  x[i] == dDNA[1]: 
                PF[1,j] = Pf[j]
    #            PB[1,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[2]: 
                PF[2,j] = Pf[j]
    #            PB[2,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[3]: 
                PF[3,j] = Pf[j]
    #            PB[3,j] = max(Pb[j], Pb_[j])     
            if j == 3:
                nfn = np.array(list([0]*len(dDNA)))
                nfd = np.array(list([0]*len(dDNA)))
#            nbn_ = np.array(list([0]*len(LDNA)))
#            nbd_ = np.array(list([0]*len(LDNA)))
#            nbn = np.array(list([0]*len(LDNA)))
#            nbd = np.array(list([0]*len(LDNA))) 
        P0 = np.dot([1,1,1,1], PF)
        if l == 0:
            t.append(P0)
            childsN.append(4)
            l += 1
            for j in range(len(dDNA)):
                parents.append(dDNA[j])
                t.append(PF[j,:])
                childsN.append(4)        
                for j in range(len(dDNA)):
                    if i - 1 ==  i - Lf +1 :
                        nfn[j] = Seq[i-1].count(dDNA[j])
                    else:
                        nfn[j] = Seq[i-Lf+1:i-1+1].count(dDNA[j])
                    if nfn[j] == 0:
                        nfn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
                    if i  ==  i - Lf:
                        nfd[j] = Seq[i].count(dDNA[j])
                    else:
                        nfd[j] = Seq[i-Lf+1:i+1].count(dDNA[j])
                    if nfd[j] == 0:
                        nfd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)+0.00001 
                    if nfd[j] < 0.001:
                        if nfn[j] == 0:
                            Pf[j] = 0
                        else:
                            Pf[j] = 1
                    else:
                        Pf[j] = nfn[j]/nfd[j]  
                
                    KL = sum(Pf*np.log(2, P0/Pf))   
                    if KL <= c:
                        Pf[j] = prune(Pf[j], j, Lf, Seq, i, dDNA, l, c, 'f')           
                    if x[i] == dDNA[0]:
                        PF[0,j] = Pf[j]
                #            PB[0,j] = max(Pb[j], Pb_[j])               
                    if  x[i] == dDNA[1]: 
                        PF[1,j] = Pf[j]
                #            PB[1,j] = max(Pb[j], Pb_[j])
                    if  x[i] == dDNA[2]: 
                        PF[2,j] = Pf[j]
                #            PB[2,j] = max(Pb[j], Pb_[j])
                    if  x[i] == dDNA[3]: 
                        PF[3,j] = Pf[j]
                #            PB[3,j] = max(Pb[j], Pb_[j])     
                    if j == 3:
                        nfn = np.array(list([0]*len(dDNA)))
                        nfd = np.array(list([0]*len(dDNA)))
                temp = np.dot([1,1,1,1], PF)
                for j in range(len(dDNA)):
                    parents.append(dDNA[temp!= 0])
                    t.append(PF[0:,temp[temp!= 0]])
                l += 1

        l += 1
        
  
        
#            if max(0, k -1) == max(0, k - Lf):
#                nfn[j] = Seq[max(0, k-1)].count(LDNA[j]) + nfn[j]
#                nfn_[j] = len(x) - nfn[j] + nfn_[j]
#            else:
#                Min = min(max(0, k-1), max(0, k-Lf))
#                Max = max(max(0, k-1), max(0, k-Lf))
#                nfn[j] = Seq[Min:Max].count(LDNA(j)) + nfn[j]
#                nfn_[j] = len(x) - nfn[j] + nfn_[j]
#            if max(0, k) == max(0, k - Lf):
#                nfd[j] = Seq[max(0, k)].count(LDNA[j]) + nfd[j]
#                nfd_[j] = len(x)-nfd[j] + nfd_[j]
#            else:
#                Min = min(max(0, k), max(0, k-Lf))
#                Max = max(max(0, k), max(0, k-Lf))
#                nfd[j] = Seq[Min:Max].count(LDNA(j)) + nfd[j]
#                nfd_[j] = len(x)-nfd[j] + nfd_[j]
#            if max(0, k -1) == max(0, k - Lb):
#                nbn[j] = Seq[max(0, k-1)].count(LDNA[j]) + nbn[j]
#                nbn_[j] = len(x)-nbn[j] + nbn_[j]
#            else:
#                Min = min(max(0, k-1), max(0, k-Lb))
#                Max = max(max(0, k-1), max(0, k-Lb))
#                nbn[j] = Seq[Min:Max].count(LDNA(j)) + nbn[j]
#                nbn_[j] = len(x)-nbn[j] + nbn_[j]
#            if max(0, k) == max(0, k - Lb):
#                nbd[j] = Seq[max(0, k)].count(LDNA[j]) + nbd[j]
#                nbd_[j] = len(x)-nbd[j] + nbd_[j]
#            else:
#                Min = min(max(0, k), max(0, k-Lb))
#                Max = max(max(0, k), max(0, k-Lb))
#                nbd[j] = Seq[Min :Max].count(LDNA(j)) + nbd[j] 
#                nbd_[j] = len(x)-nbd[j] + nbd_[j]            
           
#
#    Seq[i] = max(s)
#    for j in range(len(LDNA)):
#        if max(0, i - 1) == max(0, i - Lb):
#            nbn[j] = Seq[max(0, i-1)].count(LDNA[j])
#            if nbn[j] == 0:
#                nbn[j] = Seq[max(0, i-1)].count(LDNA[j])
#            nbn_[j] = 1-nbn[j]
#        else:
#            Min = min(max(0, i-1), max(0, i-Lb))
#            Max = max(max(0, i-1), max(0, i-Lb))
#            nbn[j] = Seq[Min:Max].count(LDNA(j))
#            nbn_[j] = Lb-1-nbn[j]
#        if max(0, i) == max(0, i - Lb):
#            nbd[j] = Seq[max(0, i)].count(LDNA[j])
#            nbd_[j] = 1-nbd[j]
#        else:
#            Min = min(max(0, i), max(0, i-Lb))
#            Max = max(max(0, i), max(0, i-Lb))
#            nbd[j] = Seq[Min :Max].count(LDNA(j))
#            nbd_[j] = Lb-nbd[j]             
#       
#        Pf[j] = nfn[j]/(nfd[j]+0.0001)
#        Pf_[j] = nfn_[j]/(nfd_[j]+0.0001)
#            
#        if Pf[j]> Pf_[j]:
#            Lp0[j] = 1          
#        if  x[i] == LDNA[0]:
#            PF[0,j] = max(Pf[j], Pf_[j])
##            PB[0,j] = max(Pb[j], Pb_[j])               
#        if  x[i] == LDNA[1]: 
#            PF[1,j] = max(Pf[j], Pf_[j])
##            PB[1,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[2]: 
#            PF[2,j] = max(Pf[j], Pf_[j])
##            PB[2,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[3]: 
#            PF[3,j] = max(Pf[j], Pf_[j])
##            PB[3,j] = max(Pb[j], Pb_[j])     
#        if j == 3:
#            nfn_ = np.array(list([0]*len(LDNA)))
#            nfd_ = np.array(list([0]*len(LDNA)))
#            nfn = np.array(list([0]*len(LDNA)))
#            nfd = np.array(list([0]*len(LDNA)))
##            nbn_ = np.array(list([0]*len(LDNA)))
##            nbd_ = np.array(list([0]*len(LDNA)))
##            nbn = np.array(list([0]*len(LDNA)))
##            nbd = np.array(list([0]*len(LDNA))) 
#            Lp0 = np.array(list([0]*len(LDNA)))    
#            
#            
#        Pb[j] = nbn[j]/(nbd[j]+0.0001)
#        Pb_[j] = nbn_[j]/(nbd_[j]+0.0001)
#            
#        if Pf[j]> Pf_[j] & k == i:
#            Lp0[j] = 1          
#        if  x[i] == LDNA[0]:
#            PF[0,j] = max(Pf[j], Pf_[j])
#            PB[0,j] = max(Pb[j], Pb_[j])               
#        if  x[i] == LDNA[1]: 
#            PF[1,j] = max(Pf[j], Pf_[j])
#            PB[1,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[2]: 
#            PF[2,j] = max(Pf[j], Pf_[j])
#            PB[2,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[3]: 
#            PF[3,j] = max(Pf[j], Pf_[j])
#            PB[3,j] = max(Pb[j], Pb_[j])  
#    
#        if j == 3:
#            nfn_ = np.array(list([0]*len(LDNA)))
#            nfd_ = np.array(list([0]*len(LDNA)))
#            nfn = np.array(list([0]*len(LDNA)))
#            nfd = np.array(list([0]*len(LDNA)))
#            nbn_ = np.array(list([0]*len(LDNA)))
#            nbd_ = np.array(list([0]*len(LDNA)))
#            nbn = np.array(list([0]*len(LDNA)))
#            nbd = np.array(list([0]*len(LDNA))) 
#            Lp0 = np.array(list([0]*len(LDNA)))
        
        
        
        
    
      
                