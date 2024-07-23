import matplotlib.pyplot as plt
import numpy as np
import levitacao

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s', nome="teste") 

zm = 1e-3*np.linspace(0,10.5, 200)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))


SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
#SO.criaRefletor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3])
SO.criaBola(1e-3, 10, [0,0,10.1e-3])
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = np.pi)


P = SO.calculaMedP2(coord, 4)

fig = plt.figure(dpi=300)

plt.title("Pressão em 2ª ordem")


plt.plot(1e3*zm, P, marker='', linestyle ='-')

plt.xlabel( 'posição (mm)' )
plt.ylabel( 'pressao (Pa)')
plt.grid()     

plt.show()