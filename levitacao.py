import numpy as np
import math


class SimuladorOndas:
    _prefix = {'': 1,  'c': 1e-2,  'm': 1e-3,  'u': 1e-6,  'n': 1e-9, 'k': 1e3, 'M': 1e-6, 'G':1e-9 }
    _pSI = {'L':1, 'T':1, 'M':1e3}
    
    
    def __init__(self,rhoSI = 1.2, c0SI = 343, f0SI=40e3, uL='mm', uT='ms', uM='g', nome = ""):
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
        
        self.bolas =[]
        self.nomeBolas =[]
        
        self.nome = nome
        
        if self.nome !="":
          f = open(self.nome + '_Simulacao.txt', "w")
          f.write("rhoSI = 1.2, c0SI = 343, f0SI=40e3, uL="+uL+", uT="+uT+", uM="+uM+"\n")
          f.write("Emissores\n")
          f.write("a, da, N, P, U0, fase, nome\n")
          f.write("Refletores\n")
          f.write("a, da, N, P, nome\n")
          f.write("Bolas\n")
          f.write("r, n_pontos_eq, P, nome\n")
          f.close()
        
    def criaRefletor(self, a, da, n, P, nome=""):
        if self.nome !="":
            f = open(self.nome + '_Simulacao.txt', "r")
            ntexto = ""
            for i in range(0,5):
                ntexto = ntexto + f.readline()
    
            for i in self.emissores :
                ntexto = ntexto + f.readline()
            
            for i in self.refletores :
                ntexto = ntexto + f.readline()
                
            ntexto = ntexto +str(a)+","+str(da)+","+str(n)+","+str(P) + ","+str(nome) + "\n"
            
            linha = f.readline()
            while linha:
                ntexto = ntexto + linha
                linha = f.readline()
            f.close()
            
            f = open(self.nome + '_Simulacao.txt', "w")
            f.write(ntexto)
            f.close()
        
        self.refletores.append(self.Refletor(self, a, da, n, P))
        self.nomeRefletores.append(nome)
        
    def criaEmissor(self, a, da, n, P, U0, fase=0,nome=""):
        if self.nome !="":
            f = open(self.nome + '_Simulacao.txt', "r")
            ntexto = ""
            for i in range(0,3):
                ntexto = ntexto + f.readline()
    
            for i in self.emissores :
                ntexto = ntexto + f.readline()
    
            ntexto = ntexto +str(a)+","+str(da)+","+str(n)+","+str(P)+","+ str(U0) + ","+str(nome) + "\n"
            linha = f.readline()
            while linha:
                ntexto = ntexto + linha
                linha = f.readline()
            f.close()
            
            f = open(self.nome + '_Simulacao.txt', "w")
            f.write(ntexto)
            f.close()
        
        self.emissores.append(self.Emissor(self, a, da, n, P, U0, fase))
        self.nomeEmissores.append(nome)
        
        self.refletoresEmissores.append(self.Refletor(self, a, da, n, P))
        
        
    def criaBola(self, r, n, P0, nome = ""):
        if self.nome !="":
            f = open(self.nome + '_Simulacao.txt', "r")
            ntexto = ""
            for i in range(0,7):
                ntexto = ntexto + f.readline()
    
            for i in self.emissores :
                ntexto = ntexto + f.readline()
            
            for i in self.refletores :
                ntexto = ntexto + f.readline()
            
            for i in self.bolas :
                ntexto = ntexto + f.readline()
                
            ntexto = ntexto +str(r)+","+str(n)+","+str(P0) + ","+str(nome) + "\n"
            
            linha = f.readline()
            while linha:
                ntexto = ntexto + linha
                linha = f.readline()
            f.close()
            
            f = open(self.nome + '_Simulacao.txt', "w")
            f.write(ntexto)
            f.close()        
        
        self.bolas.append(self.Bola(self, r, n, P0))
        self.nomeBolas.append(nome)
        
    def calculaPar(self, Pts, Nref, Pressao = True, Deslocamento = True):
        refEff = np.concatenate((self.refletoresEmissores, self.refletores, self.bolas))#isso garante que a lista começa sempre pelos emissores fazendo papel de refletores; também junta as bolas como relfetores
        Ptotal = [] 
        Dtotal = []
        #T-->M
        if Pressao:
            P = np.zeros(len(Pts))
            for em in self.emissores:
               P = np.add(P, em.pressao(Pts))     
            Ptotal.append(P)
        if Deslocamento:
            D = np.zeros((len(Pts),3))
            for em in self.emissores:
               D = np.add(D, em.deslocamento(Pts))     
            Dtotal.append(D)
        
        if(Nref != 0):
            A = []
            indA = np.zeros((len(refEff), 2), dtype=np.integer)
            for i, ref in enumerate(refEff):
                A.append(ref.superficie()) 
                if i!=0:
                    indA[i,0] = indA[i-1, 1]+1
                    indA[i,1] = indA[i,0]+len(ref.superficie())-1
                else:
                    indA[i,1] = len(ref.superficie())-1
            A = np.concatenate(A) 
            
            #T-->R
            P=np.zeros(len(A))
            for T, ind in zip(self.emissores, indA):   
                a1,a2,a3 =np.split(A,[ind[0],1+ind[1]], axis=0)
                p1 = T.pressao(a1)
                p2 = np.zeros(len(a2))
                p3 = T.pressao(a3)
                P = P + np.concatenate((p1, p2, p3))
            Pint = P #guarda o valor da pressão intermediária na superficie dos refletores
            
            #---------------------
            
            #R-->M                
            if Pressao:
                P = np.zeros(len(Pts))
                for R, ind in zip(refEff, indA):
                   P = P + R.pressao(Pts,Pint[ind[0]:ind[1]+1])     
                Ptotal.append(P)
            
            if Deslocamento:
                D = np.zeros((len(Pts),3))
                for R, ind in zip(refEff, indA):
                   D = D + R.deslocamento(Pts,Pint[ind[0]:ind[1]+1])  
                Dtotal.append(D)
                
                  
            for N in range (0,Nref-1):
                #R-->R
                
                P=np.zeros(len(A))
                for R, ind in zip(refEff, indA):   
                    a1,a2,a3 =np.split(A,[ind[0],1+ind[1]], axis=0)
                    p1 = R.pressao(a1, Pint[ind[0]:1+ind[1]])
                    p2 = np.zeros(len(a2))
                    p3 = R.pressao(a3, Pint[ind[0]:1+ind[1]])
                    P = P + np.concatenate((p1, p2, p3))
                Pint = P #guarda o valor da pressão intermediária na superficie dos refletores
               
                #R-->M                
                if Pressao:
                    P = np.zeros(len(Pts))
                    for R, ind in zip(refEff, indA):
                       P = P + R.pressao(Pts,Pint[ind[0]:ind[1]+1])     
                    Ptotal.append(P)
                
                if Deslocamento:
                    D = np.zeros((len(Pts),3))
                    for R, ind in zip(refEff, indA):
                       D = D + R.deslocamento(Pts,Pint[ind[0]:ind[1]+1])  
                    Dtotal.append(D)
                
        if self.nome != "":
            ca = "x , y, z"
            if Pressao:
                for i in range (0, len(Ptotal)):
                    ca = ca +  ", Pnref " + str(i)
            if Deslocamento:
                for i in range (0, len(Dtotal)):
                    ca = ca +  ", Dxnref " + str(i) + ", Dynref " + str(i) + ", Dznref " + str(i)
                
            Dprint = np.transpose(Dtotal, axes = (0,2,1))
            Dprint = np.concatenate(Dprint)
            dados = np.array([*np.transpose(Pts), *Ptotal, *Dprint])
            dados = np.transpose(dados)
            
            np.savetxt(self.nome + '_pressao.csv',dados, header=ca, delimiter=',')
            print("dados salvos em .csv")

        return Ptotal, Dtotal
    
    def calculaPeD(self, Pts, Nref):
        refEff = np.concatenate((self.refletoresEmissores, self.refletores, self.bolas))#isso garante que a lista começa sempre pelos emissores fazendo papel de refletores; também junta as bolas como relfetores
        
        PDtotal =[]
        #T-->M
       
        PD = np.zeros((len(Pts),4))
        for em in self.emissores:
           PD = np.add(PD, em.PeD(Pts))     
        PDtotal.append(PD)
        
        if(Nref != 0):
            A = []
            indA = np.zeros((len(refEff), 2), dtype=np.integer)
            for i, ref in enumerate(refEff):
                A.append(ref.superficie()) 
                if i!=0:
                    indA[i,0] = indA[i-1, 1]+1
                    indA[i,1] = indA[i,0]+len(ref.superficie())-1
                else:
                    indA[i,1] = len(ref.superficie())-1
            A = np.concatenate(A)
                    
            #T-->R

            P=np.zeros(len(A))
            for T, ind in zip(self.emissores, indA):   
                a1,a2,a3 =np.split(A,[ind[0],1+ind[1]], axis=0)
                p1 = T.pressao(a1)
                p2 = np.zeros(len(a2))
                p3 = T.pressao(a3)
                P = P + np.concatenate((p1, p2, p3))
            Pint = P #guarda o valor da pressão intermediária na superficie dos refletores
            
            #---------------------
            
            #R-->M
    
            PD = np.zeros((len(Pts),4))
            for R, ind in zip(refEff, indA):
               PD = PD + R.PeD(Pts, Pint[ind[0]:ind[1]+1])
            PDtotal.append(PD)
            
            for N in range (0,Nref-1):
                #R-->R
                
                P=np.zeros(len(A))
                for R, ind in zip(refEff, indA):   
                    a1,a2,a3 =np.split(A,[ind[0],1+ind[1]], axis=0)
                    p1 = R.pressao(a1, Pint[ind[0]:1+ind[1]])
                    p2 = np.zeros(len(a2))
                    p3 = R.pressao(a3, Pint[ind[0]:1+ind[1]])
                    P = P + np.concatenate((p1, p2, p3))
                Pint = P #guarda o valor da pressão intermediária na superficie dos refletores
                
                #R-->M
                
                PD = np.zeros((len(Pts),4))
                for R, ind in zip(refEff, indA):
                   PD = PD + R.PeD(Pts, Pint[ind[0]:ind[1]+1])
                PDtotal.append(PD)
                
        if self.nome != "":
            ca = "x , y, z"

            for i in range (0, Nref+1):
                ca = ca +  ", Pnref " + str(i)+  ", Dxnref " + str(i) + ", Dynref " + str(i) + ", Dznref " + str(i)
                
            PDprint = np.transpose(PDtotal, axes = (0,2,1))
            PDprint = np.concatenate(PDprint)
            dados = np.array([*np.transpose(Pts), *PDprint])
            dados = np.transpose(dados)
            
            np.savetxt(self.nome + '_pressao.csv',dados, header=ca, delimiter=',')
            print("dados salvos em .csv")

        return PDtotal
    
    def calculaPar2Bola(self, PtsMed, PtsBo, Nref, Raio):
        refEff = np.concatenate((self.refletoresEmissores, self.refletores))#isso garante que a lista começa sempre pelos emissores fazendo papel de refletores
            
        PintTotal=[]
        
        if(Nref != 0):
            A = []
            indA = np.zeros((len(refEff), 2), dtype=np.integer)
            for i, ref in enumerate(refEff):
                A.append(ref.superficie()) 
                if i!=0:
                    indA[i,0] = indA[i-1, 1]+1
                    indA[i,1] = indA[i,0]+len(ref.superficie())-1
                else:
                    indA[i,1] = len(ref.superficie())-1
            A = np.concatenate(A)
                    
            #T-->R

            P=np.zeros(len(A))
            for T, ind in zip(self.emissores, indA):   
                a1,a2,a3 =np.split(A,[ind[0],1+ind[1]], axis=0)
                p1 = T.pressao(a1)
                p2 = np.zeros(len(a2))
                p3 = T.pressao(a3)
                P = P + np.concatenate((p1, p2, p3))
            Pint = P #guarda o valor da pressão intermediária na superficie dos refletores
            PintTotal.append(P)
            
            
            for N in range (0,Nref-1):
                #R-->R
                
                P=np.zeros(len(A))
                for R, ind in zip(refEff, indA):   
                    a1,a2,a3 =np.split(A,[ind[0],1+ind[1]], axis=0)
                    p1 = R.pressao(a1, Pint[ind[0]:1+ind[1]])
                    p2 = np.zeros(len(a2))
                    p3 = R.pressao(a3, Pint[ind[0]:1+ind[1]])
                    P = P + np.concatenate((p1, p2, p3))
                Pint = P #guarda o valor da pressão intermediária na superficie dos refletores
                PintTotal.append(P)
            
            PrefBol = np.sum(np.delete(PintTotal, len(PintTotal)-1, axis=0), axis = 0)
                                       

        PrefMed = np.sum(PintTotal, axis = 0)
        
           
        
        PDr = np.zeros((len(PtsMed),4))
        Pr = np.zeros(len(PtsBo))
        if (Nref != 0):
            for R, ind in zip(refEff, indA):
               PDr = PDr + R.PeD(PtsMed, PrefMed[ind[0]:ind[1]+1])
              
            for R, ind in zip(refEff, indA):
               Pr = Pr + R.pressao(PtsBo, PrefBol[ind[0]:ind[1]+1])
        
        
        
        PDt = np.zeros((len(PtsMed),4))
        for em in self.emissores:
           PDt = PDt + em.PeD(PtsMed)     
        Pt = np.zeros(len(PtsBo))
        for em in self.emissores:
           Pt = Pt + em.pressao(PtsBo) 
           
        PeD0med = PDt+PDr
        PBol = Pt+Pr
      
        PeD1Med = np.zeros((len(PtsBo), len(PtsMed), 4), dtype = np.complexfloating)
        
        k = 2*math.pi*self.f/self.c0
        A = 4*math.pi*(Raio**2)
        multV = (complex(0,-1)/((self.c0/self.f)*k*self.rho*self.c0))
        multP = complex(0,1)/(self.c0/self.f)
        for i, med in enumerate(PtsMed):
            
            for j, bo in enumerate(PtsBo):
                rlinha = np.linalg.norm(med-bo)
                rlinhadir = (med-bo)/rlinha
                P = PBol[j]*(1/rlinha)*(math.e**(complex(0,-1)*k*rlinha))*A
                V = rlinhadir*((complex(0,1)/rlinha)-k)*P         
                
                P = P*multP
                V = V*multV
                
                PeD1Med[j,i, :] = PeD0med[i, :] + (P, *V)
                
                
                
        P = np.absolute(PeD1Med[:,:,0])
        D = np.absolute(np.delete(PeD1Med, 0, axis = 2))
        
        P2med = (P**2)/2 
        D2med = np.sum(D**2, axis=2)/2
        
        
        Gorkov =(P2med/(2*self.rho*(self.c0**2))) -(D2med*self.rho*3/4)
        
        return Gorkov
             
    
    def calculaMedP2(self, Pts, Nref):
        PDref = self.calculaPeD(Pts, Nref)
        
        PD = np.absolute(np.sum(PDref, axis=0))
        P = PD[:,0]
        D = np.delete(PD, 0, 1)
        P2med = (P**2)/2 
        D2med = np.sum(D**2, axis=1)/2
      
        medP2 =(P2med/(2*self.rho*(self.c0**2))) -(D2med*self.rho/2)
        
        if self.nome != "":
            ca = "x, y, z, medP2"
           
            dados = np.array([*np.transpose(Pts), medP2])
            dados = np.transpose(dados)
            
            np.savetxt(self.nome + '_pressao2ordem.csv',dados, header=ca, delimiter=',')
            print("dados salvos em .csv")

        return medP2
    
    def calculaGorkov(self, Pts, Nref):
        PDref = self.calculaPeD(Pts, Nref)
        
        PD = np.absolute(np.sum(PDref, axis=0))
        P = PD[:,0]
        D = np.delete(PD, 0, 1)
        P2med = (P**2)/2 
        D2med = np.sum(D**2, axis=1)/2
      
        Gorkov =(P2med/(2*self.rho*(self.c0**2))) -(D2med*self.rho*3/4)
        
        if self.nome != "":
            ca = "x, y, z, Gorkov"
           
            dados = np.array([*np.transpose(Pts), Gorkov])
            dados = np.transpose(dados)
            
            np.savetxt(self.nome + '_Gorkov.csv',dados, header=ca, delimiter=',')
            print("dados salvos em .csv")

        return Gorkov
    
    def reCalculaMedP2(self, P, D):  #função que calcula a pressão em segunda ordem, porém já recebendo os valores de pressão e deslocamento em 1 ordem      

        P2med = (np.array(P)**2)/2   
        D2med = np.sum(np.array(D)**2, axis=1)/2
            
        medP2 =(P2med/(2*self.rho*(self.c0**2))) -(D2med*self.rho/2)

        return medP2
            
        
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
        
        def PeD(self, Pts0):
            Pts = np.array(Pts0)
            PD =[]
            for Pt in Pts:
                if (np.dot((Pt-self.Pem[0]), self.N)>0): #checa se o ponto  está na frente do emissor
                    Pi=0
                    Vi = 0
                    for dr in  self.Pem:
                        rlinha = np.linalg.norm(Pt-dr)
                        rlinhadir = (Pt-dr)/rlinha
                        dP = (1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf
                        Pi = Pi+dP
                        dV = rlinhadir*((complex(0,1)/rlinha)-self.k)*dP                
                        Vi = Vi+dV
                    Pi = (complex(0,1)*self.rho*self.c0*self.U/self.Lamb)*Pi
                    Vi = (complex(0,-1)*self.U/(self.Lamb*self.k))*Vi
                    
                    PD.append((Pi, *Vi))
                else:
                    PD.append([complex(0,0)*4])
            return PD
            
    class Bola():
        
        def __init__(self, outer, r, n, P0):
            self.f = outer.f
            self.c0 = outer.c0
            self.rho = outer.rho
           
            self.w = 2*math.pi*self.f
            self.k = self.w/self.c0
            self.Lamb = self.c0/self.f
            
            self.P=P0            
            if n > 1:
                self.Pbo = []
                self.A = []
                theta  = np.linspace(math.pi/n,math.pi, n-1)
                for t in theta:
                    phi = np.linspace(0,math.pi*2, math.ceil(n*math.sin(t)))
                    dA = (r**2)*(math.pi/n)*(math.pi*2/math.ceil(n*math.sin(t)))*math.sin(t)
                    for ph in phi:
                        dP = [r*math.sin(t)*math.cos(ph)+self.P[0],r*math.sin(t)*math.sin(ph)+self.P[1], r*math.cos(t)+self.P[2]]
                        self.Pbo.append(np.array(dP))
                        self.A.append(dA)
                #correção da área:
                A=0
                for dA in self.A:
                    A = A+ dA
                k = (4*math.pi*(r**2))/A
                self.A = np.multiply(k, self.A)
            else:
                self.Pbo = [self.P]
                self.A = [(4*math.pi*(r**2))]
        

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
        
        def PeD(self,Pts0, Pinc):
            Pts = np.array(Pts0)
            if len(self.Pbo)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            PD=[]

            for Pt in Pts:
                Pi=0
                Vi=0
                for dr, da, Pin in  zip(self.Pbo, self.A, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    rlinhadir = (Pt-dr)/rlinha
                    dP = Pin*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*da
                    Pi = Pi+dP
                    dV = rlinhadir*((complex(0,1)/rlinha)-self.k)*dP
                    Vi = Vi+dV
                Pi = Pi*complex(0,1)/self.Lamb
      
                Vi = Vi*(complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))
                PD.append((Pi, *Vi)) 
            return PD

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
        
        def PeD(self,Pts0, Pinc):
            Pts = np.array(Pts0)
            if len(self.Prf)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            PD=[]
            
            for Pt in Pts:
                Pi=0
                Vi=0
                for dr, Pin in  zip(self.Prf, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    rlinhadir = (Pt-dr)/rlinha
                    dP = Pin*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf
                    Pi = Pi+dP
                    dV = rlinhadir*((complex(0,1)/rlinha)-self.k)*dP                   
                    Vi = Vi+dV
                Pi = Pi*complex(0,1)/self.Lamb
                  
                Vi = Vi*(complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))
                PD.append((Pi, *Vi))
            return PD
        
        

              


