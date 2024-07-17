import numpy as np
import math
import matplotlib.pyplot as plt

class SimuladorOndas:
    _prefix = {'': 1,  'c': 1e-2,  'm': 1e-3,  'u': 1e-6,  'n': 1e-9, 'k': 1e3, 'M': 1e-6, 'G':1e-9 }
    _pSI = {'L':1, 'T':1, 'M':1e3}
    
    
    def __init__(self,rhoSI = 1.2, c0SI = 343, f0SI=40e3, uL='mm', uT='ms', uM='g'):
        if(len(uL)==1): self.pL=""
        else: self.pL = uL[0]
        if(len(uM)==1): self.pM=""
        else: self.pM = uM[0]
        if(len(uT)==1): self.pT="" 
        else: self.pT = uT[0]
 
        self.rho = rhoSI*(self._pSI["M"]/self._prefix[self.pM])*(self._pSI["L"]/self._prefix[self.pL])**(-3)
        self.c0 = c0SI*(self._pSI["M"]/self._prefix[self.pM])*(self._pSI["T"]/self._prefix[self.pT])**(-1)
        self.f = f0SI*(self._pSI["T"]/self._prefix[self.pT])**(-1)
        
        self.refletores = []
        self.nomeRefletores = []
        
        self.emissores = []
        self.refletoresEmissores = []
        self.nomeEmissores = []
        
        
        
    def criaRefletor(self, a, da, n, P, nome=""):
        self.refletores.append(self.Refletor(self, a, da, n, P))
        self.nomeRefletores.append(nome)
        
    def criaEmissor(self, a, da, n, P, U0, fase=0,nome=""):
        self.emissores.append(self.Emissor(self, a, da, n, P, U0, fase))
        self.nomeEmissores.append(nome)
        
        self.refletoresEmissores.append(self.Refletor(self, a, da, n, P))
    
    def criaBola(self, r, n, P0, nome = ""):
        self.refletores.append(self.Bola(self, r, n, P0))
        self.nomeRefletores.append(nome)
        
    def calculaP(self, Pts, Nref):
        refEff = np.concatenate((self.refletoresEmissores, self.refletores))#isso garante que a lista começa sempre pelos emissores fazendo papel de refletores
        Ptotal = [] 
        P = np.zeros(len(Pts))
        for em in self.emissores:
           P = np.add(P, em.pressao(Pts))     
        Ptotal.append(P)
        if(Nref != 0):
            A = []
            for ref in refEff:
               A.append(ref.superficie()) 
               Pint = []#guarda o valor da pressão intermediária na superficie dos refletores
            #T-->R
            P = np.zeros(A.shape())
            for T, i in zip(self.emissores, range(0, len(self.emissores))):
                PT =[]
                for Ai, j in zip(A, range(0,len(A))):
                    if(j!=i):
                        PT.append(T.pressao(Ai))
                    else:
                        PT.append(np.zeros(Ai.shape()))
                P = np.add(P, PT)
                
            #R-->R

            for R, i in zip(refEff, range(0, len(refEff))):
                PT =[]
                for Ai, j in zip(A, range(0,len(A))):
                    if(j!=i):
                        PT.append(R.pressao(Ai))
                    else:
                        PT.append(np.zeros(Ai.shape()))
                P = np.add(P, PT)
        return Ptotal
        
    class Emissor():
        
        def __init__(self, outer, a, da, n, P, U0, fase=0):
            self.f = outer.f
            self.c0 = outer.c0
            self.rho = outer.rho
            
            self.U = U0*(math.e**(complex(0,1)*fase))
            
            self.w = 2*math.pi*self.f
            self.k = self.w/self.c0
            self.Lamb = self.c0/self.f
            
            
            
            if (math.floor(2*a/da)>1):
                xi = np.linspace(-1,1,math.floor(2*a/da), True)
                yi = xi
            else:
                xi = np.array([0])
                yi=xi
            #cria os vetores geradores do plano do emissor
            v1 = np.cross(n,[1,0,0])
            if np.linalg.norm(v1) !=0:
                v1 = (a/np.linalg.norm(v1))*v1
                v2 = np.cross(n, v1)/np.linalg.norm(n)
            else:
                v1 = np.array([0,a,0])
                v2 = np.array([0,0,a])
            
            self.N =n/np.linalg.norm(n)#guarda o valor da normal como vetor unitário
            self.Pem = []
            #cria todos os pontos do emissor (eles são uma combinação linear de v1 e v2)
            #e seu raio é menor que a
            for x in xi:
                for y in yi:
                    rd = (v1*x) + (v2*y)
                    if (np.linalg.norm(rd) <= a):
                        self.Pem.append(P+rd)
            
            if(len(self.Pem)>1):
                self.dAEf = (math.pi*(a**2))/len(self.Pem)
            else:
                self.Pem=[P]
                self.dAEf = math.pi*(a**2)
            
        def pressao(self,Pts0):
            Pts = np.array(Pts0)
            P=[]
            for Pt in Pts:
                if (np.dot((Pt-self.Pem[0]), self.N)>0): #checa se o ponto  está na frente do emissor
                    Pi=0
                    for dr in  self.Pem:
                        rlinha = np.linalg.norm(Pt-dr)
                        dP = (1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf
                        Pi = Pi+dP
                    Pi = (complex(0,1)*self.rho*self.c0*self.U/self.Lamb)*Pi
                    P.append(Pi)         
                else:
                    P.append(complex(0,0))
            return P
        
        def deslocamento(self, Pts0):
            Pts = np.array(Pts0)
            V=[]
            for Pt in Pts:
                dist = np.dot((Pt-self.Pem[0]), self.N)
                if (dist>0): #checa se o ponto  a ser calculado a velocidade está na frente do emissor
                    Vi=0
                    for dr in  self.Pem:
                        rlinha = np.linalg.norm(Pt-dr)
                        rlinhadir = (Pt-dr)/rlinha
                        dV = ((complex(0,1)/rlinha)-self.k)*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf                    
                        dV = dV*rlinhadir
                        Vi = Vi+dV
                    Vi = (complex(0,-1)*self.U/(self.Lamb*self.k))*Vi
                    V.append(Vi) 
                elif dist ==0:
                    V.append(self.U*self.N)
                else:
                    V.append([0,0,0])
            return V
        
    class Bola():
        
        def __init__(self, outer, r, n, P0):
            self.f = outer.f
            self.c0 = outer.c0
            self.rho = outer.rho
           
            self.w = 2*math.pi*self.f
            self.k = self.w/self.c0
            self.Lamb = self.c0/self.f
            
            self.P=P0            
            self.Pbo = []
            self.A = []
            theta  = np.linspace(math.pi/n,math.pi, n-1)
            for t in theta:
                phi = np.linspace(0,math.pi*2, math.ceil(n*math.sin(t)))
                dA = (r**2)*(math.pi/n)*(math.pi*2/math.ceil(n*math.sin(t)))*math.sin(t)
                for ph in phi:
                    dP = [r*math.sin(t)*math.cos(ph)+self.P[0],r*math.sin(t)*math.sin(ph)+self.P[1], r*math.cos(t)+self.P[2]]
                    self.Pbo.append(dP)
                    self.A.append(dA)
            #correção da área:
            A=0
            for dA in self.A:
                A = dA
            k = (4*math.pi*r**2)/A
            for dA in self.A:
                dA=k*dA

        def superficie(self):
            return self.Pbo
         

        
        def pressao(self,Pts0, Pinc):
            Pts = np.array(Pts0)
            if len(self.Pbo)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            P=[]
            for Pt in Pts:
                Pi=0
                for dr, da, Pin in  zip(self.Pbo, self.A, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    dP = Pin*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*da
                    Pi = Pi+dP
                Pi = Pi*complex(0,1)/self.Lamb
                P.append(Pi)         
            return P
        
        def deslocamento(self, Pts0,Pinc):
            Pts = np.array(Pts0)
            if len(self.Pbo)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            V=[]
            for Pt in Pts:
                Vi=0
                for dr, da, Pin in  zip(self.Pbo,self.A, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    rlinhadir = (Pt-dr)/rlinha
                    dV = ((complex(0,1)/rlinha)-self.k)*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*da                  
                    dV = Pin*dV*rlinhadir
                    Vi = Vi+dV
                Vi = Vi*(complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))
                V.append(Vi) 
                
            return V

    class Refletor():
        
        def __init__(self, outer, a, da, n, P):
            self.f = outer.f
            self.c0 = outer.c0
            self.rho = outer.rho
            
            self.w = 2*math.pi*self.f
            self.k = self.w/self.c0
            self.Lamb = self.c0/self.f
            
            
            
            if (math.floor(2*a/da)>1):
                xi = np.linspace(-1,1,math.floor(2*a/da), True)
                yi = xi
            else:
                xi = np.array([0])
                yi=xi
            #cria os vetores geradores do plano do emissor
            v1 = np.cross(n,[1,0,0])
            if np.linalg.norm(v1) !=0:
                v1 = (a/np.linalg.norm(v1))*v1
                v2 = np.cross(n, v1)/np.linalg.norm(n)
            else:
                v1 = np.array([0,a,0])
                v2 = np.array([0,0,a])
            
            self.N =n/np.linalg.norm(n)#guarda o valor da normal como vetor unitário
            self.Prf = []
            #cria todos os pontos do emissor (eles são uma combinação linear de v1 e v2)
            #e seu raio é menor que a
            for x in xi:
                for y in yi:
                    rd = (v1*x) + (v2*y)
                    if (np.linalg.norm(rd) <= a):
                        self.Prf.append(P+rd)
            
            if(len(self.Prf)>1):
                self.dAEf = (math.pi*(a**2))/len(self.Prf)
            else:
                self.Prf=[P]
                self.dAEf = math.pi*(a**2)
            
        def superficie(self):
            return self.Prf
            
        def pressao(self,Pts0, Pinc):
            Pts = np.array(Pts0)
            if len(self.Prf)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            P=[]
            for Pt in Pts:
                Pi=0
                for dr, Pin in  zip(self.Prf, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    dP = Pin*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf
                    Pi = Pi+dP
                Pi = Pi*complex(0,1)/self.Lamb
                P.append(Pi)         
            return P
        
        def deslocamento(self, Pts0,Pinc):
            Pts = np.array(Pts0)
            if len(self.Prf)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            V=[]
            for Pt in Pts:
                Vi=0
                for dr, Pin in  zip(self.Prf, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    rlinhadir = (Pt-dr)/rlinha
                    dV = ((complex(0,1)/rlinha)-self.k)*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf                    
                    dV = Pin*dV*rlinhadir
                    Vi = Vi+dV
                Vi = Vi*(complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))
                V.append(Vi) 
                
            return V

              

SO = SimuladorOndas(uL='m', uM='kg', uT='s') 

zm = 1e-3*np.linspace(1,20, 200)
coord=[]
for zi in zm:
    coord.append(np.array([0,0,zi]))

SO.criaEmissor(7e-3, 1e-3, [0,0,1], [0,0,-1e-3], 1)
SO.criaEmissor(7e-3, 1e-3, [0,0,-1], [0,0,21e-3], 1, fase = math.pi)

P = SO.calculaP(coord, 0)
P0 = P[0]

fig = plt.figure(dpi=300)

plt.title("pressao absoluta eixo z")
plt.plot(1e3*zm, np.absolute(P0), marker='', linestyle ='-', label="Emissor")
plt.xlabel( 'posição (mm)' )
plt.ylabel( 'pressao (Pa)')
plt.grid()
plt.legend()      

plt.show()

