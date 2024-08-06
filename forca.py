import matplotlib.pyplot as plt
import numpy as np
import levitacao

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s', nome="teste") 

zm = 1e-3*np.linspace(0,20, 200)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))


SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = np.pi)


F = SO.calculaForca(coord, 0)
G = SO.calculaGorkov(coord, 0)


fig = plt.figure(dpi=300)

plt.title("Força no eixo z")

plt.subplot(2,1,1)
plt.plot(1e3*zm, F[:,2], marker='', linestyle ='-', label = "Força")
plt.subplot(2,1,2)
plt.plot(1e3*zm, G, marker='', linestyle ='-', label = "Gor'kov")

plt.xlabel( 'posição (mm)' )
plt.ylabel( 'pressao (Pa)')
plt.grid()
plt.legend()      

plt.show()