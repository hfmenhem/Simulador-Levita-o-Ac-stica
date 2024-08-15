import matplotlib.pyplot as plt
import numpy as np
import levitacao
from scipy.optimize import curve_fit

import timeit


def seno(x, A, B, fase, C):
    return A*np.sin(B*x+fase)+C

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s') 

z0 = 12e-3

zm = 1e-3*np.linspace(4,8.5, 100)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))
    
bm = 1e-3*np.linspace(-1.5,1.5, 30)
coordB=[]
for bi in bm:
    coordB.append(np.array([0,0,bi]))
    
zm0 = 1e-3*np.linspace(0,20, 100)
coord0=[]
for zi in zm0:
    coord0.append(np.array([0,0,zi]))


SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = 0)


#SO.criaRefletor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3])
SO.criaBola(.5e-3, 1, [0,0,z0])

#--medir tempo--
start = timeit.default_timer()

F0, G0 = SO.calculaForca(coord0, 6,CalGorcov=True, Bolas=False)

F, G = SO.calculaPar2Bola(coord, coordB, 6, CalGorcov=True)

stop = timeit.default_timer()
execution_time = stop - start

print("Programa demorou "+str(execution_time)+" s para calcular") # It returns time in seconds



parametros =[]
#bounds0 = [[0, 0, 0], [np.inf, 2*np.pi/20e-3, 2*np.pi]]
for Fi in F:
    Amp = (np.max(Fi[:,2])-np.min(Fi[:,2]))/2
 
    B = abs(np.pi/(zm[np.argmax(Fi[:,2])]-zm[np.argmin(Fi[:,2])]))
    bounds0 = [[0, 0, 0, -np.inf], [2*Amp, 2*B, 2*np.pi, np.inf]]
    Fpar, Fpar_cov = curve_fit(seno, zm,Fi[:,2], bounds=bounds0)
    parametros.append(Fpar)
    
equilibrio = []
grad =[]
for Pari in parametros:
    eq = (np.pi + np.arcsin(Pari[3]/Pari[0])-Pari[2])/Pari[1]
    gr = -Pari[1]*np.sqrt((Pari[0]**2)+(Pari[3]**2))
    equilibrio.append(eq)
    grad.append(gr)

equilibrio = np.array(equilibrio)
grad = np.array(grad)


plt.figure(dpi =300)
plt.plot(1e3*zm0, 1e-3*F0[:,2], marker='', label = "sem bola", linestyle ='-')

plt.xlabel( 'Posição (mm)' )
plt.ylabel( 'Força (kN/m^3)')
plt.grid()
plt.title("Sem bola")
plt.show()

plt.figure(dpi=300)
plt.subplot(2,1,1)
plt.plot(1e3*bm, 1e3*equilibrio, marker='.', linestyle ='-')
plt.ylabel( 'Posição de equilibrio (mm)')

plt.subplot(2,1,2)
plt.plot(1e3*bm, grad, marker='.', linestyle ='-')
plt.xlabel( 'Posição da Bola 2 (mm)' )
plt.ylabel( 'Gradiente Força (N/m^3/m)')
plt.show()


