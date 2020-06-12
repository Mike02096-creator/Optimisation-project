# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:13:25 2019

"""

import numpy as np

g = 9.81 #Intensité de la pesanteur (m/s²)
ro = 1.28 #Masse volumique de l'air (kg/m^3)

def lectureVitesse(fichier) : #Lecture du fichier de vitesse
    import csv
    with open(fichier, newline='') as csvfile :
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        reader = csv.reader(csvfile)
        v = []
        for row in reader :
            v += row
    vCycle = np.array([0.]*len(v))
    for i in range(len(v)) :
        vCycle[i] = float(v[i])/3.6
    return vCycle

class Scooter :
    def __init__(self,M,Cr,Ap,Cd,v) :
        self.M = M #Masse volumique de l'air (kg/m^3)
        self.Cr = Cr #Coefficient de roulement
        self.Ap = Ap #Aire frontale (m²)
        self.Cd = Cd #Coefficient de traînée
        self.v = v #Vitesse du scooter (m/s)
    def aScooter(self) : #Accélération du scooter (m/s²)
        a = []
        for i in range(len(self.v)) :
            a.append(self.v[i]-self.v[i-1])
        return np.array(a)
    def fRoulement(self) : #Résistance au roulement (N)
        return self.Cr*self.M*g
    def fAerodynamique(self) : #Résistance aérodynamique (N)
        return (ro*self.Ap*self.Cd*self.v**2)/2
    def fInertie(self) : #Inertie (N)
        return self.M*self.aScooter()
    def fTraction(self) : #Force de traction (N)
        return self.fRoulement()+self.fAerodynamique()+self.fInertie()
    def pMeca(self) : #Puissance mécanique (W)
        return self.fTraction()*self.v

class MachineElectrique :
    def __init__(self,rend,Scooter) :
        self.rend = rend #Rendement constant
        self.s = Scooter
    def pElec(self) : #Puissance électrique (W)
        return self.rend*self.s.pMeca()

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

class PileACombustible :
    def __init__(self,MachineElectrique,k) :
        self.m = MachineElectrique
        self.k = k #Coefficient de proportionnalité entre la puissance de la pile à combustile et la puissance électrique
    def pPAC(self) : #Puissance de la pile à combustible
        return self.k*self.m.pElec()
    def pChim(self) : #Puissance chimique (W)
        return c0+c1*self.pPAC()+c2*self.pPAC()**2+c3*self.pPAC()**3+c4*self.pPAC()**4
    def dmH2(self) : #Débit masssique d'H2
        return (MH2*self.pChim())/PCS
    def vPAC(self) : #Tension aux bornes de la pile (V)
        if (self.m.pElec() >= 0) and (self.m.pElec() <= pLim) :
            return a1*self.m.pElec()+b1
        else :
            return a2*self.m.pElec()+b2

class Batterie :
    def __init__(self,MachineElectrique,E0,Rint,Q0) :
        self.m = MachineElectrique
        self.E0 = E0 #Force électromotrice (V)
        self.Rint = Rint #Résistance interne (Ohm)
        self.Q0 = Q0 #Capacité (A.h)
    def iBatt(self) : #Courant (A)
        I = 0*self.m.pElec()
        delta = self.E0**2-4*self.Rint*self.m.pElec()
        for i in range(len(self.m.pElec())) :
            if delta[i] >= 0 :
                I[i] = (self.E0-np.sqrt(delta[i]))/(2*self.Rint)
            else :
                print("Erreur calcul")
                I[i] = 0
        return I
    def Q(self) : #Charge (A.h)
        q = []
        for i in range(len(self.iBatt())) :
            q.append(np.trapz(self.iBatt()[:i]))
        return self.Q0-np.array(q)
    def SOC(self) : #État de charge (%)
         return (self.Q()/self.Q0)*100
