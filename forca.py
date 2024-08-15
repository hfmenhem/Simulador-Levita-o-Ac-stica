import matplotlib.pyplot as plt
import numpy as np
import levitacao

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s', nome="teste") 

zm = 1e-3*np.linspace(0,20, 200)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))


SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = 0)


F, G = SO.calculaForca(coord, 2, CalGorcov=True)

fig = plt.figure(dpi=300)

plt.subplots_adjust(hspace=0.4)
plt.subplot(2,1,1)
plt.plot(1e3*zm, F[:,2], marker='', linestyle ='-', label = "Força")
plt.ylabel( 'Força (N/m^3) (Pa)')
plt.grid()
plt.title("Força")
plt.subplot(2,1,2)


plt.plot(1e3*zm, G, marker='', linestyle ='-', label = "Gor'kov")
plt.xlabel( 'posição (mm)' )
plt.ylabel( 'Potencial (J/m^3)')
plt.title("Potencial de Gor'kov")
plt.grid()

plt.show()