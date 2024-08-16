import matplotlib.pyplot as plt
import numpy as np
import levitacao
from scipy.optimize import curve_fit

import timeit

D = 20e-3
disc = 1e-3
R = 7e-3


SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s') 

SO.criaEmissor(R, disc, [0,0,1], [0,0,0], 1)
SO.criaEmissor(R, disc, [0,0,-1], [0,0,D], 1, fase = 0)
#SO.criaRefletor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3])

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


#Sabendo onde são os pontos de equilibrio, é possível posicionar a bola
Zb0 = zeros[1]
Zan = zeros[2]
Zans = np.linspace(Zan-(zeros[1]-zeros[0])/3, Zan+(zeros[1]-zeros[0])/3, 100)
coordAns = np.zeros([len(Zans), 3])
coordAns[:,2] = Zans

SO.criaBola(1e-3, 20, [0,0,Zb0])

Zvar = np.linspace((zeros[1]-zeros[0])/8, (zeros[1]-zeros[0])/-8, 20)
coordVar = np.zeros([len(Zvar), 3])
coordVar[:,2] = Zvar

F, G = SO.calculaPar2Bola(coordAns, coordVar, 6, CalGorcov=True)
#acha pts de equilibrio deslocados

zerosNeq =[]
gradsNeq=[]
for Fi in F[:,:,2]:
    zNeq = Zans
    FNeq = Fi

    zteste = zNeq[:-1]
    multTeste = FNeq[1:]*FNeq[:-1]
    difTeste = FNeq[1:]-FNeq[:-1]
    condNeq = np.logical_and(multTeste<0, difTeste<0 )
    intNeq = np.transpose(((zteste[:-1])[condNeq[:-1]],  (zteste[1:])[condNeq[:-1]]))
    FintNeq = np.transpose(((FNeq[:-2])[condNeq[:-1]],  (FNeq[1:-1])[condNeq[:-1]]))
    
    zeroNeq = (intNeq[:,0]*FintNeq[:,1]-intNeq[:,1]*FintNeq[:,0])/(FintNeq[:,1]-FintNeq[:,0])
    gradNeq = (FintNeq[:,1]-FintNeq[:,0])/(intNeq[:,1]-intNeq[:,0])
    if len(zNeq ==1):
        zerosNeq.append(*zeroNeq)
        gradsNeq.append(*gradNeq)
    else:
        print("há mais de 1 ponto de equilíbrio")
    
zerosNeq = np.array(zerosNeq)
gradsNeq = np.array(gradsNeq) 
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


#Gráfico da perturbação
plt.figure(dpi=300)

plt.subplot(2,1,1)

plt.title("Com 1 bola variando")
plt.plot(1e3*(Zvar+Zb0), 1e3*zerosNeq, marker='.', linestyle ='-')
plt.plot([1e3*(Zvar[0]+Zb0), 1e3*(Zvar[-1]+Zb0)], [1e3*Zan, 1e3*Zan], marker='', linestyle ='-')

plt.ylabel( 'Posição de equilibrio \n (mm)')
plt.grid()

plt.subplots_adjust(hspace=0.4)

plt.subplot(2,1,2)

plt.plot(1e3*(Zvar+Zb0), 1e-6*gradsNeq, marker='.', linestyle ='-')
#plt.plot([1e3*(bm[1]+z0), 1e3*(bm[-1]+z0)], [1e-6*grad0, 1e-6*grad0], marker='', linestyle ='-')

plt.xlabel( 'Posição da Bola 2 (mm)' )
plt.ylabel( 'Gradiente Força \n (kN/m^3/nm)')
plt.grid()
plt.show()