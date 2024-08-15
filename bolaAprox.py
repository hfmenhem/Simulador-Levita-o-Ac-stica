import matplotlib.pyplot as plt
import numpy as np
import levitacao


SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s') 
z0 = 12e-3
zm = 1e-3*np.linspace(6.5,8.5, 100)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))
    
bm = 1e-3*np.linspace(-1.5,1.5, 10)
coordB=[]
for bi in bm:
    coordB.append(np.array([0,0,bi]))


SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = 0)


#SO.criaRefletor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3])
SO.criaBola(.5e-3, 1, [0,0,z0])

F0, G0 = SO.calculaForca(coord, 2,CalGorcov=True, Bolas=False)
F, G = SO.calculaPar2Bola(coord, coordB, 2, CalGorcov=True)


plt.figure(dpi =300)
plt.subplots_adjust(hspace=0.4)
plt.subplot(2,1,1)
for i, Fi in enumerate(F[:,:,2]):
    plt.plot(1e3*zm, 1e-3*Fi, marker='', label = str(np.round((bm[i]+z0)*1e3, decimals=1)) +"mm", linestyle ='-')
plt.plot(1e3*zm, 1e-3*F0[:,2], marker='', label = "sem bola", linestyle ='-')

plt.xlabel( 'posição (mm)' )
plt.ylabel( 'Força (kN/m^3)')
plt.grid()
plt.legend(loc='upper left',bbox_to_anchor=(1, 1))  

plt.subplot(2,1,2)
for i, Gi in enumerate(G):
    plt.plot(1e3*zm, Gi, marker='', linestyle ='-')
plt.plot(1e3*zm, G0, marker='', label = "sem bola", linestyle ='-')
plt.ylabel( 'Potencial (J/m^3)')
plt.grid()


plt.show()