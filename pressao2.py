import matplotlib.pyplot as plt
import numpy as np
import levitacao

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s') 

zm = 1e-3*np.linspace(6,9, 200)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))
    
zm2 = 1e-3*np.linspace(11, 14, 5)
coord2=[]
for zi in zm2:
    coord2.append(np.array([0,0,zi]))


SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
#SO.criaRefletor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3])
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = 0)

#SO.criaBola(.5e-3, 1, [0,0,12e-3])

A = SO.calculaPar2Bola(coord, coord2, 4, 0.5e-3)
G = SO.calculaGorkov(coord, 4)
#PD = np.sum(SO.calculaPeD(coord, 4), axis=0)

#P = SO.calculaMedP2(coord, 4)

fig = plt.figure(dpi=300)
plt.title("Gorkov")
for i in range(0, len(A)):    
    plt.plot(1e3*zm, A[i,:], marker='', linestyle ='-', label = str(round(1e3*coord2[i][2], 1))  +"mm")
    
plt.plot(1e3*zm, G, marker='', linestyle ='-', label = "sem bola")

plt.xlabel( 'posição (mm)' )
plt.ylabel( 'Gorkov/volume  (J/m^3)')
plt.grid()    
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))   

plt.show()