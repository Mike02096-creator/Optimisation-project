# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 13:47:56 2020

@author: Daniel
"""

import numpy as np
import matplotlib.pyplot as plt
from projet import lectureVitesse, Scooter, MachineElectrique, PileACombustible, Batterie

g = 9.81 #Intensité de la pesanteur (m/s²)
ro = 1.28 #Masse volumique de l'air (kg/m^3)

#Coefficients puissance chimique
c0 = 84.50
c1 = 2.185
c2 = 0
c3 = -9.715e-8
c4 = 1.722e-9

PCS = 283e3 #Pouvoir calorifique supérieur de l'H2 (J/mol)
MH2 = 2.016e-3 #Masse molaire de l'H2 (kg/mol)

#Coefficients tension
a1 = -0.0528
b1 = 40.7721
a2 = -0.0175
b2 = 37.5178

pLim = 100 #Puissance limite (W)

vS = lectureVitesse('Cycle_WLTC.csv') #Vitesse du scooter (m/s)
S = Scooter(110,0.01,0.7,0.75,vS)
M = MachineElectrique(0.8,S)
k1 = 0 #Coefficient du cas i
PAC1 = PileACombustible(M,k1)
B = Batterie(M,12,2*10**(-3),20)

u1 = PAC1.pPAC() #Commande
w = M.pElec() #Perturbation
x = B.SOC() #Variable d'état
t0 = 0 #Temps initial (s)
tf = len(vS)-1 #Temps final (s)
t = np.arange(t0,tf+1,1) #Intervalle de temps

def JH2_1(u) : #Coût
    return np.trapz((MH2*(c0+c1*u+c2*u**2+c3*u**3+c4*u**4)/PCS)[t0:tf])

lambda0 = 0 #Borne supérieure de l'intervalle de lambda
lambda1 = -0.9728353140916808 #Borne inférieure de l'intervalle de lambda
lbd = np.linspace(lambda1,lambda0,len(vS)) #Intervalle de lambda
eps = 10**(-5) #Précision

def f(u) : #Contrainte
    return -(100*B.E0-np.sqrt(B.E0**2-4*B.Rint*(w-u)))/(2*B.Q0*B.Rint)

f1 = f(u1)

def H(u) : #Hamiltonien
    mat = []
    for i in f(u) :
        ligne = []
        for j in lbd :
            ligne.append(i*j)
        mat.append(ligne)
    return PAC1.dmH2()+np.array(mat)

def dichotomie(x0,xf,eps,n) :
    x=[0] * 5    
    y=[0] * 5      
    x[0]=x0
    x[4]=xf   
    i=0    
    while abs(x[4]-x[0]) > eps :     
        x[2]=(x[0]+x[4])/2
        x[1]=(x[0]+x[2])/2
        x[3]=(x[2]+x[4])/2
        for j in range (5) :            
            y[j]=H(x[j])[:,n][0]
        if y[0] < y[1] and y[1] < y[2] and y[2] < y[3] and y[3] < y[4] :
            x[4] = x[1]
        elif y[0] > y[1] and y[1] < y[2] and y[2] < y[3] and y[3] < y[4] :
            x[4] = x[2]
        elif y[0] > y[1] and y[1] > y[2] and y[2] < y[3] and y[3] < y[4] :
            x[0] = x[1]
            x[4] = x[3]
        elif y[0] > y[1] and y[1] > y[2] and y[2] > y[3] and y[3] < y[4] :
            x[0] = x[2]
        elif y[0] > y[1] and y[1] > y[2] and y[2] > y[3] and y[3] > y[4] :
            x[0] = x[3]
        i += 1
    xmin = x[0]
    for k in x :
        if H(k)[:,n][0] < H(xmin)[:,n][0] :
            xmin = k
    return xmin

def PMP(eps,u,lbd,t0,tf) : #Algorithme de Pontryaguine
    i = -1
    d = 0
    while i < len(lbd)-1 and x[tf] != x[t0] :
        i += 1
        d = dichotomie(u[0],u[len(u)-1],eps,i)
        tf -= 1
    return lbd[i],d

plt.figure()
plt.plot(t,S.v)
plt.figure()
plt.plot(t,S.aScooter())

PMP1 = PMP(eps,u1,lbd,t0,tf)
print("lambda* =",PMP1[0],"; u* =",PMP[1])
print("La fonction coût minimisée vaut",JH2_1(np.array([PMP1[1]]*(len(t)+1))))