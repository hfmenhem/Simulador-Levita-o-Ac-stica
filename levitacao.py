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
    
    def calculaPeD(self, Pts, Nref, forca = False, Bolas = True, Pintermediaria = False):
        if Bolas:
            refEff = np.concatenate((self.refletoresEmissores, self.refletores, self.bolas))#isso garante que a lista começa sempre pelos emissores fazendo papel de refletores; também junta as bolas como relfetores
        else:
            refEff = np.concatenate((self.refletoresEmissores, self.refletores))
        PDtotal =[]
        #T-->M
       
        if forca:
           Par = np.zeros((len(Pts),13))
           for em in self.emissores:
              Par = np.add(Par, em.ParForca(Pts))     
           PDtotal.append(Par) 
        else:
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
            if forca:
               Par = np.zeros((len(Pts),13))
               for R, ind in zip(refEff, indA):
                  Par = Par + R.ParForca(Pts, Pint[ind[0]:ind[1]+1])
               PDtotal.append(Par)
            else:
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
                
                if forca:
                   Par = np.zeros((len(Pts),13))
                   for R, ind in zip(refEff, indA):
                      Par = Par + R.ParForca(Pts, Pint[ind[0]:ind[1]+1])
                   PDtotal.append(Par)
                else:
                    PD = np.zeros((len(Pts),4))
                    for R, ind in zip(refEff, indA):
                       PD = PD + R.PeD(Pts, Pint[ind[0]:ind[1]+1])
                    PDtotal.append(PD)
                
        if self.nome != "":
            ca = "x , y, z"
            
            if forca:
                for i in range (0, Nref+1):
                    ca = ca +  ", Pnref " + str(i)+  ", Dxnref " + str(i) + ", Dynref " + str(i) + ", Dznref " + str(i) + ", dUx/dxNref " + str(i) + ", dUy/dxNref " + str(i)+ ", dUz/dxNref " + str(i)+ ", dUx/dyNref " + str(i) + ", dUy/dyNref " + str(i)+ ", dUz/dyNref " + str(i) + ", dUx/dzNref " + str(i) + ", dUy/dzNref " + str(i)+ ", dUz/dzNref " + str(i)
            else:
                for i in range (0, Nref+1):
                    ca = ca +  ", Pnref " + str(i)+  ", Dxnref " + str(i) + ", Dynref " + str(i) + ", Dznref " + str(i)
                
            PDprint = np.transpose(PDtotal, axes = (0,2,1))
            PDprint = np.concatenate(PDprint)
            dados = np.array([*np.transpose(Pts), *PDprint])
            dados = np.transpose(dados)
            
            np.savetxt(self.nome + '_pressao.csv',dados, header=ca, delimiter=',')
            print("dados salvos em .csv")
        
        if Pintermediaria:
            return PDtotal, Pint
        else:
            return PDtotal
    
    def calculaPar2Bola(self, PtsMed, PtsBo, Nref, CalGorcov=False):
        PtsBo = np.transpose(np.array([PtsBo]), (1,0,2))
        PtsMed = np.array(PtsMed)

        A0 =np.array( [self.bolas[0].superficie()])       
        Ab = A0 + PtsBo
        
        #calculo dos valores base para os pontos de medida e bolas
        Parl ,Pint = self.calculaPeD(np.concatenate([PtsMed, np.reshape(Ab, (-1,3))]), Nref-1, forca = True, Bolas = False, Pintermediaria = True)
        
        Par = np.sum(Parl, axis=0)
        PM = Par[:int(np.size(np.array(PtsMed))/3), :]
        Pbl = Par[int(np.size(np.array(PtsMed))/3):, :]
        
        Pb = np.reshape(Pbl, [*np.shape(Ab[:,:,0]), 13])

    
        #calculo da última reflexão para os pontos de medida:
        refEff = np.concatenate((self.refletoresEmissores, self.refletores))
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
        
        Pmf = np.zeros((len(PtsMed),13))
        for R, ind in zip(refEff, indA):
           Pmf = Pmf + R.ParForca(PtsMed, Pint[ind[0]:ind[1]+1])
        PM = PM + Pmf
        
        #calculo da pressão vinda da bola:
            
        #B-->M
        
        Pfb =[]
        for pinc, d in zip(Pb[:,:,0], PtsBo):
           Pfb.append(self.bolas[0].ParForca(PtsMed, pinc))
        Pfb = np.array(Pfb)
        
        ParFinal = Pfb + np.array([PM])
        
        #calcula força e Gor'kov
        
        k = 2*math.pi*self.f/self.c0
        

        P = ParFinal[:,:,[0]]
        D = ParFinal[:,:, 1:4]
        GDx = ParFinal[:,:, 4:7]
        GDy = ParFinal[:,:, 7:10]
        GDz = ParFinal[:,:, 10:]
        
        p1 = np.imag(np.multiply(np.conjugate(P), D))
        p2 = np.real((np.conjugate(D[:,:,[0]])*GDx) + (np.conjugate(D[:,:,[1]])*GDy) + (np.conjugate(D[:,:,[2]])*GDz))
        
        
        Forca = -1*((p1*k/(2*self.c0)) -(p2*self.rho*3/4))

        if CalGorcov:
            P2med = (np.absolute(P[:,:,0])**2)/2 
            D2med = np.sum(np.absolute(D)**2, axis=2)/2
          
            Gorkov =(P2med/(2*self.rho*(self.c0**2))) -(D2med*self.rho*3/4)
            return [Forca, Gorkov]
        else:
            return Forca
        return
             
    
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
    
    def calculaForca(self, Pts, Nref, CalGorcov = False, Bolas = True):
        Parref = self.calculaPeD(Pts, Nref,forca = True, Bolas=Bolas)
        k = 2*math.pi*self.f/self.c0
        
        Par = np.sum(Parref, axis=0)

        P = Par[:,[0]]
        D = Par[:, 1:4]
        GDx = Par[:, 4:7]
        GDy = Par[:, 7:10]
        GDz = Par[:, 10:]
        
        p1 = np.imag(np.multiply(np.conjugate(P), D))
        p2 = np.real((np.conjugate(D[:,[0]])*GDx) + (np.conjugate(D[:,[1]])*GDy) + (np.conjugate(D[:,[2]])*GDz))
        
        Forca = -1*((p1*k/(2*self.c0)) -(p2*self.rho*3/4))
        
        if CalGorcov:
            P2med = (np.absolute(P[:,0])**2)/2 
            D2med = np.sum(np.absolute(D)**2, axis=1)/2
          
            Gorkov =(P2med/(2*self.rho*(self.c0**2))) -(D2med*self.rho*3/4)
        
        if self.nome != "":
            if CalGorcov:
                ca = "x, y, z, Fx, Fy, Fz, Gor'kov"
                dados = np.array([*np.transpose(Pts), *np.transpose(Forca), Gorkov])
            else:
                ca = "x, y, z, Fx, Fy, Fz"
                dados = np.array([*np.transpose(Pts), *np.transpose(Forca)])
            
            dados = np.transpose(dados)
            
            np.savetxt(self.nome + '_Forca.csv',dados, header=ca, delimiter=',')
            print("dados salvos em .csv")

        if CalGorcov:
            return [Forca, Gorkov]
        else:
            return Forca
    
        
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
        
        def ParForca(self, Pts0):
            Pts = np.array(Pts0)
            Par =[]
            for Pt in Pts:
                if (np.dot((Pt-self.Pem[0]), self.N)>0): #checa se o ponto  está na frente do emissor
                    Pi=0
                    Vi = 0
                    Vdxi = 0
                    Vdyi = 0
                    Vdzi = 0
                    for dr in  self.Pem:
                        rvec = Pt-dr
                        rlinha = np.linalg.norm(rvec)
                        rlinhadir = (Pt-dr)/rlinha
                        
                        C1 = (1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf
                        C2 = ((complex(0,1)/rlinha)-self.k)
                        C3 = C1/rlinha                      
                        C4 = rlinhadir*C3*((complex(0,-3)/(rlinha**2))+(3*self.k/rlinha)+(complex(0,1)*(self.k**2)))
                        C5 = C3*C2
                        
                        dP = C1
                        dV = rlinhadir*C1*C2
                        dVdxi = [C5,0,0] + C4*rvec[0]
                        dVdyi = [0,C5,0] + C4*rvec[1]
                        dVdzi = [0,0,C5] + C4*rvec[2]
                                                
                        Pi = Pi+dP               
                        Vi = Vi+dV
                        Vdxi = Vdxi + dVdxi
                        Vdyi = Vdyi + dVdyi
                        Vdzi = Vdzi + dVdzi
                        
                    Pi = (complex(0,1)*self.rho*self.c0*self.U/self.Lamb)*Pi
                    Vi = (complex(0,-1)*self.U/(self.Lamb*self.k))*Vi
                    Vdxi = (complex(0,-1)*self.U/(self.Lamb*self.k))*Vdxi
                    Vdyi = (complex(0,-1)*self.U/(self.Lamb*self.k))*Vdyi
                    Vdzi = (complex(0,-1)*self.U/(self.Lamb*self.k))*Vdzi
                    
                    Par.append((Pi, *Vi, *Vdxi, *Vdyi, *Vdzi))
                else:
                    Par.append([complex(0,0)*13])
            return Par
            
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
                        dP = [r*math.sin(t)*math.cos(ph) + self.P[0],r*math.sin(t)*math.sin(ph) + self.P[1], r*math.cos(t)+self.P[2]]
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
            

        def superficie(self, Dcentro = 0):
            return np.add(self.Pbo,  Dcentro)
        
        def centro(self):
            return self.P
        
        def pressao(self,Pts0, Pinc, Dcentro = 0):
            PboC = np.add(self.Pbo,  Dcentro)
            Pts = np.array(Pts0)
            if len(PboC)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            P=[]
            for Pt in Pts:
                Pi=0
                for dr, da, Pin in  zip(PboC, self.A, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    dP = Pin*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*da
                    Pi = Pi+dP
                Pi = Pi*complex(0,1)/self.Lamb
                P.append(Pi)         
            return P
        
        def deslocamento(self, Pts0,Pinc, Dcentro = 0):
            PboC = np.add(self.Pbo,  Dcentro)
            Pts = np.array(Pts0)
            if len(PboC)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            V=[]
            for Pt in Pts:
                Vi=0
                for dr, da, Pin in  zip(PboC,self.A, Pinc):
                    rlinha = np.linalg.norm(Pt-dr)
                    rlinhadir = (Pt-dr)/rlinha
                    dV = ((complex(0,1)/rlinha)-self.k)*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*da                  
                    dV = Pin*dV*rlinhadir
                    Vi = Vi+dV
                Vi = Vi*(complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))
                V.append(Vi) 
                
            return V
        
        def PeD(self,Pts0, Pinc, Dcentro = 0):
            PboC = np.add(self.Pbo,  Dcentro)
            Pts = np.array(Pts0)
            if len(PboC)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            PD=[]

            for Pt in Pts:
                Pi=0
                Vi=0
                for dr, da, Pin in  zip(PboC, self.A, Pinc):
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
        
        def ParForca(self,Pts0, Pinc,Dcentro = 0):
            PboC = np.add(self.Pbo,  Dcentro)
            Pts = np.array(Pts0)
            if len(PboC)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            Par =[]
            for Pt in Pts:
                Pi=0
                Vi = 0
                Vdxi = 0
                Vdyi = 0
                Vdzi = 0
                for dr, da, Pin in  zip(PboC, self.A, Pinc):
                    rvec = Pt-dr
                    rlinha = np.linalg.norm(rvec)
                    rlinhadir = (Pt-dr)/rlinha
                    
                    C1 = Pin*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*da
                    C2 = ((complex(0,1)/rlinha)-self.k)
                    C3 = C1/rlinha                      
                    C4 = rlinhadir*C3*((complex(0,-3)/(rlinha**2))+(3*self.k/rlinha)+(complex(0,1)*(self.k**2)))
                    C5 = C3*C2
                    
                    dP = C1
                    dV = rlinhadir*C1*C2
                    dVdxi = [C5,0,0] + C4*rvec[0]
                    dVdyi = [0,C5,0] + C4*rvec[1]
                    dVdzi = [0,0,C5] + C4*rvec[2]
                                            
                    Pi = Pi+dP               
                    Vi = Vi+dV
                    Vdxi = Vdxi + dVdxi
                    Vdyi = Vdyi + dVdyi
                    Vdzi = Vdzi + dVdzi
                    
                Pi = Pi*complex(0,1)/self.Lamb
                Vi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vi
                Vdxi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vdxi
                Vdyi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vdyi
                Vdzi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vdzi
                
                Par.append((Pi, *Vi, *Vdxi, *Vdyi, *Vdzi))
                
            return Par

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
            return np.array(self.Prf)
            
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
        
        def ParForca(self,Pts0, Pinc):
            Pts = np.array(Pts0)
            if len(self.Prf)!=len(Pinc):
                raise Exception("tamanho de Pem deve ser igual ao de Pinc")
            Par =[]
            for Pt in Pts:
         
                Pi=0
                Vi = 0
                Vdxi = 0
                Vdyi = 0
                Vdzi = 0
                
                for dr, Pin in  zip(self.Prf, Pinc):
                    rvec = Pt-dr
                    rlinha = np.linalg.norm(rvec)
                    rlinhadir = (Pt-dr)/rlinha
                    
                    C1 = Pin*(1/rlinha)*(math.e**(complex(0,-1)*self.k*rlinha))*self.dAEf
                    C2 = ((complex(0,1)/rlinha)-self.k)
                    C3 = C1/rlinha                      
                    C4 = rlinhadir*C3*((complex(0,-3)/(rlinha**2))+(3*self.k/rlinha)+(complex(0,1)*(self.k**2)))
                    C5 = C3*C2
                    
                    dP = C1
                    dV = rlinhadir*C1*C2
                    dVdxi = [C5,0,0] + C4*rvec[0]
                    dVdyi = [0,C5,0] + C4*rvec[1]
                    dVdzi = [0,0,C5] + C4*rvec[2]
                                            
                    Pi = Pi+dP               
                    Vi = Vi+dV
                    Vdxi = Vdxi + dVdxi
                    Vdyi = Vdyi + dVdyi
                    Vdzi = Vdzi + dVdzi
                    
                Pi = Pi*complex(0,1)/self.Lamb
                Vi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vi
                Vdxi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vdxi
                Vdyi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vdyi
                Vdzi = (complex(0,-1)/(self.Lamb*self.k*self.rho*self.c0))*Vdzi
                
                Par.append((Pi, *Vi, *Vdxi, *Vdyi, *Vdzi))
                
            return Par
        
        

              


