import matplotlib.pyplot as plt
import numpy as np
import levitacao

SO = levitacao.SimuladorOndas(uL='m', uM='kg', uT='s', nome="teste") 

zm = 1e-3*np.linspace(0,11, 200)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))


SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
#SO.criaRefletor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3])
SO.criaBola(1e-3, 10, [0,0,12e-3])
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = np.pi)


#P, D = SO.calculaPar(coord, 4, Pressao=False, Deslocamento=True)
PD = SO.calculaPeD(coord, 4)
Pl = []
for p, i in zip(PD, range(0, len(PD))):
    if i == 0:
        Pl.append(np.transpose(p)[0])
    else:
        Pl.append(Pl[i-1]+np.transpose(p)[0])
      

Dlz = []
for d, i in zip(PD, range(0, len(PD))):
    if i == 0:
        Dlz.append(np.transpose(d)[3])
    else:
        Dlz.append(Dlz[i-1]+np.transpose(d)[3])


fig = plt.figure(dpi=300)

plt.title("velocidade absoluta eixo z")

for i in range(0, len(Pl)):
    plt.plot(1e3*zm, np.absolute(Pl[i]), marker='', linestyle ='-', label= str(i) + " nref")

plt.xlabel( 'posição (mm)' )
plt.ylabel( 'pressao (Pa)')
plt.grid()
plt.legend()      

plt.show()