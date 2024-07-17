import matplotlib.pyplot as plt
import numpy as np
import levitacao

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s') 

zm = 1e-3*np.linspace(1,20, 200)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))

SO.criaEmissor(7e-3, .5e-3, [0,0,1], [0,0,-1e-3], 1)
SO.criaRefletor(7e-3, .5e-3, [0,0,-1], [0,0,21e-3])
#SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = math.pi)

P = SO.calculaP(coord, 4)

Pl = []
for p, i in zip(P, range(0, len(P))):
    if i == 0:
        Pl.append(p)
    else:
        Pl.append(Pl[i-1]+p)

fig = plt.figure(dpi=300)

plt.title("pressao absoluta eixo z")

for i in range(0, len(Pl)):
    plt.plot(1e3*zm, np.absolute(Pl[i]), marker='', linestyle ='-', label= str(i) + " nref")

plt.xlabel( 'posição (mm)' )
plt.ylabel( 'pressao (Pa)')
plt.grid()
plt.legend()      

plt.show()