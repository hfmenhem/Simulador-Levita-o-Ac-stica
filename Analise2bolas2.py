import matplotlib.pyplot as plt
import numpy as np
import levitacao
from scipy.optimize import curve_fit

import timeit

D = 20e-3
disc = 1e-3
R = 7e-3
Zb = 12e-3

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s') 

SO.criaEmissor(R, disc, [0,0,1], [0,0,0], 1)
SO.criaEmissor(R, disc, [0,0,-1], [0,0,D], 1, fase = 0)
#SO.criaRefletor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3])

SO.criaBola(1e-3, 20, [0,0,Zb])

zT = np.linspace(disc, D-disc, 100)
coordT = np.zeros([len(zT), 3])
coordT[:,2] = zT

F0, G0 = SO.calculaForca(coordT, 6,CalGorcov=True, Bolas=False)
F0z = F0[:,2]

#Acha pt equilibrio
zZeros = zT
FZeros = F0z
for i in range(0,2):
    zteste = zZeros[:-1]
    multTeste = FZeros[1:]*FZeros[:-1]
    difTeste = FZeros[1:]-FZeros[:-1]
    cond0 = np.logical_and(multTeste<0, difTeste<0 )
    int0 = np.transpose(((zteste[:-1])[cond0[:-1]],  (zteste[1:])[cond0[:-1]]))
    zZeros =np.ndarray.flatten(np.transpose(np.linspace(int0[:,0], int0[:,1])))  
    coordzeros = np.zeros([len(zZeros), 3])
    coordzeros[:,2] =zZeros
    FZeros, GZeros = SO.calculaForca(coordzeros, 6,CalGorcov=True, Bolas=False)
    FZeros=FZeros[:,2]
zteste = zZeros[:-1]
multTeste = FZeros[1:]*FZeros[:-1]
difTeste = FZeros[1:]-FZeros[:-1]
cond0 = np.logical_and(multTeste<0, difTeste<0 )
int0 = np.transpose(((zteste[:-1])[cond0[:-1]],  (zteste[1:])[cond0[:-1]]))
Fint0 = np.transpose(((FZeros[:-2])[cond0[:-1]],  (FZeros[1:-1])[cond0[:-1]]))

zeros = (int0[:,0]*Fint0[:,1]-int0[:,1]*Fint0[:,0])/(Fint0[:,1]-Fint0[:,0])
grad = (Fint0[:,1]-Fint0[:,0])/(int0[:,1]-int0[:,0])

#gráfico Total

plt.figure(dpi =300)
plt.plot(1e3*zT, 1e-3*F0z, marker='', label = "sem bola", linestyle ='-')
z0s = np.ndarray.flatten(int0)
f0s = np.ndarray.flatten(Fint0)
plt.plot(1e3*z0s, 1e-3*f0s, marker='.', label = "sem bola", linestyle ='')
plt.plot(1e3*zeros, np.zeros(len(zeros)), marker='.', label = "sem bola", linestyle ='')
#plt.plot(1e3*zZeros, 1e-3*FZeros, marker='.', label = "sem bola", linestyle ='')

plt.xlim(0, D*1e3)
plt.xlabel( 'Posição (mm)' )
plt.ylabel( 'Força (kN/m^3)')
plt.title("Sem bola")
plt.grid()
plt.show()