import numpy as np

#Rigidez Local de Elemento:

#Montagem da matriz constitutiva EPT de um elemento com E e v
def MontaCEPT(E,v):
    c = E/(1-v**2)
    
    Cv = np.array([[c,c*v,0],[c*v,c,0,],[0,0,c*(1-v)/2]])
    
    return Cv

def MontadN(r,s):
    #Matriz com as derivadas em relação a r (primeira linha) e s(segunda linha)
    
    dNrs = (0.25)*np.array([[s-1,1-s,1+s,-(1+s)],[r-1,-(1+r),1+r,1-r]])
    
    return dNrs

def MontaJ(dNrs,X,Y):
    J = np.zeros((2,2))
    
    for i in range(X.size):
        J[0,0] += X[i]*dNrs[0,i]
        J[0,1] += Y[i]*dNrs[0,i]
        J[1,0] += X[i]*dNrs[1,i]
        J[1,1] += Y[i]*dNrs[1,i]
    
    return J

def CorrigedN(dNrs,J):

    dNxy = (np.linalg.inv(J))@dNrs
    
    return dNxy

def MontaB(dNxy):
    B = np.zeros((3,8))
    
    c = 0
    
    for i in range(4):
        B[0,c] += dNxy[0,i]
        B[1,c+1] += dNxy[1,i]
        B[2,c] += B[1,c+1]
        B[2,c+1] += B[0,c]
        
        c += 2
    return B

#Monta Ke

def MontaKe(X,Y,E,v,te):
    # Aloca Matriz
    Ke = np.zeros((8,8))
    #Matriz Constitutiva do elemento
    
    Cv = MontaCEPT(E,v) # OK
    
    #Define os pontos de Gauss
    
    pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])
    
    #Laço pelos ponos de Gauss
    for p in range(4):
        r = pg[0,p]
        s = pg[1,p]
        
        #Calcula as derivadas da função de interpolação
        dNrs = MontadN(r,s)
        
        #Calcula a Matriz Jacobiana
        J = MontaJ(dNrs,X,Y)
        
        #Calcula determinante
        dJ = np.linalg.det(J)
        
        #Mapeia as derivadas para xy
        dNxy = CorrigedN(dNrs,J)
        
        #Monta B
        B = MontaB(dNxy)

        #Acumula na matriz de rigidez do elemento

        Ke += ((np.transpose(B))@Cv@B)*te*dJ
    return Ke

def MontaXY(e,IJ,XY):
    nos = IJ[e-1]
    
    X = np.zeros((nos.size))
    Y = np.zeros((nos.size))
    
    #Preenche os vetores
    
    for i in range(nos.size):
        #Nó
        no = IJ[e-1,i]
        #Coordenadas
        X[i] = XY[no-1,0]
        Y[i] = XY[no-1,1]
        
    return X,Y

def Montagls(e,IJ):
    #Inicializa vetor
    nos = IJ[e-1]

    gls = np.zeros((2*nos.size,1))      
    c = 0
    
    for i in range(nos.size):
        #Nó
        no = nos[i]
        #Laço pelos graus de liberdade dos nós:
        for j in range(2):
            gls[c] =  int(2*(no-1)+j)
            c+=1
    
    return gls
        
    
def RigidezGlobal(nn,ne,MAT,ESP,XY,IJ,bolha = False):
    K = np.zeros((2*nn,2*nn))
    
    #Laço nos elementos
    for e in range(ne):
        Ee = MAT[e,0]
        ve = MAT[e,1]
        te = ESP[e]
        
        #Determina os gls globais do elemento:
        g = Montagls(e,IJ)
        
        # Recupera as coordenadas do elemento:
        X,Y = MontaXY(e,IJ,XY)
        
        #Monta Ke
        if g.size == 6:
            Ke = MontaKeTriangulo(X,Y,Ee,ve,te)
        elif g.size == 8:
            if bolha:
                Ke = MontaKeBolha(X,Y,Ee,ve,te)
            else:
                Ke = MontaKe(X,Y,Ee,ve,te)
        for i in range(g.size):
                for j in range(g.size):
                    K[int(g[i]),int(g[j])] += Ke[i,j]
    return K

def Faces(face):
    #Retorna dois índices para o cálculo dos vetores locais da face
    
    if face == 1 :
        p1 = 1
        p2 = 2
    elif face == 2:
        p1 =2
        p2 =3    
    elif face == 3:
        p1 =3
        p2=4
    elif face == 4:
        p1=4
        p2=1
    return p1,p2

def MontaN(ζ,p1,p2):
    #Aloca N
    N = np.zeros((2,8))
    
    #Funções de interpolação não-nulas
    
    N1 = (1-ζ)/2
    N2 = (1+ζ)/2
    
    #Posiciona na Matriz N
    N[0,2*p1-2] = N1
    N[1,2*p1-1] = N1
    N[0,2*p2-2] = N2
    N[1,2*p2-1] = N2
    
    return N

def Facesnv(p1,p2,e,IJ,XY):
    #Recupera os nós da face
    n1 = IJ[e-1,p1-1]
    n2 = IJ[e-1,p2-1]
    
    #Recupera as coordenadas
    x1 = XY[n1-1,0]
    y1 = XY[n1-1,1]
    x2 = XY[n2-1,0]
    y2 = XY[n2-1,1]
    
    #Deltas
    deltax = x2-x1
    deltay = y2-y1
    
    L = (deltax**2 + deltay**2)**(0.5)
    #Determinante do Jacobiano
    dJ = L/2
    
    # Vetor tangente à face
    v = (1/L) * np.array([[deltax],[deltay]])
    # Vetor normal à face
    n = (1/L) * np.array([[deltay],[-deltax]])
    
    return dJ,n,v

def Ftbilinear(dJ,d,p1,p2,valor,te):
    #Pontos de gauss
    
    pg = (1/(3**(0.5))) * np.array([[-1],[1]])
    
    #Aloca o vetor
    
    F = np.zeros((8,1))
    
    #Laço nos pontos de Gauss:
    
    for i in range(2):
        ζ = pg[i]
        #Monta matriz N
        N =MontaN(ζ,p1,p2)
        
        F = F+ valor*(N.transpose()@d)*dJ*te
        
    return F

#Montagem do vetor global de forças distribuídas nas faces dos elementos e diretamente nos nós
def ForcapGlobal(nn,ESP,XY,IJ,nc,P):
    F = np.zeros((2*nn,1))
    
    for i in range(nc):
        #Recupera os dados do carregamento
        ele = P[i,0]
        face = P[i,1]
        dire = P[i,2]
        val = P[i,3]

        ## Força num ponto
        if dire == 0:
            no_global = ele
            gdl_local = face
            F[2*(no_global-1)+(gdl_local-1)] += val    
        
        #Força numa face
        else:
            #Determina os gls do elemento
            g =  Montagls(ele,IJ)
            #Espessura do elemento
            te = ESP[ele-1]
            
            if g.size == 8:
                #Ponteiros da face
                p1,p2 = Faces(face)
                
                #Vetores da face
                dJ,n,v = Facesnv(p1,p2,ele,IJ,XY)
                
                if dire == 1:
                    d = n
                elif dire == 2:
                    d = v
                    
                #Vetor de Força
                Fe = Ftbilinear(dJ,d,p1,p2,val,te)
                
            elif g.size == 6:
                p1,p2 = FacesTriangulo(face)
                #Vetores da face
                dJ,n,v = FacesnvTriangulo(p1,p2,ele,IJ,XY)
                
                if dire == 1:
                    d = n
                elif dire == 2:
                    d = v
                Fe = FtTriangulo(dJ,d,p1,p2,val,te)
            
            
            for j in range(g.size):
                    F[int(g[j])] += Fe[j]
    return F
        
#Retorna o vetor com as tensões em um ponto (r,s) de um elemento bilinear isoparamétrico de 4 nós.

def CalculaTensaoElemento(r,s,X,Y,E,v,Ue):
    #Matriz constitutiva do elemento
    Cv = MontaCEPT(E,v)
    #Calcula as derivadas das funções de interpolação
    dNrs = MontadN(r,s)
    #Calcula Jacobiana
    J = MontaJ(dNrs,X,Y)    
    #Mapeia as derivadas para xy
    dNxy = CorrigedN(dNrs,J)
    
    #Monta B
    B = MontaB(dNxy)
    
    sigmae = Cv@B@Ue
    
    return sigmae
    

def CalculaTensaoMalha(nn,ne,MAT,ESP,XY,IJ,U,bolha=False):  
    
    Sigma = np.zeros((3*ne,1))
    
    for e in range(ne):
        
        Ee = MAT[e,0]
        ve = MAT[e,1]
                
        # Recupera as coordenadas do elemento:
        X,Y = MontaXY(e,IJ,XY)
        # Determina os gdls globais do elemento
        g = Montagls(e+1,IJ)
        
        Ue = np.zeros((g.size,1))

        for i in range(g.size):
            Ue[i] = U[int(g[i])]
        
        if g.size == 8:
            if bolha:
                te = ESP[e]
                sigmae = CalculaTensaoElementoBolha(0,0,X,Y,Ee,ve,Ue,te)
            else:
                sigmae = CalculaTensaoElemento(0,0,X,Y,Ee,ve,Ue)
        elif g.size == 6:
            sigmae = CalculaTensaoElementoTriangulo(0,0,X,Y,Ee,ve,Ue)
        for j in range(3):
            Sigma[3*e+j] = sigmae[j]
        
    
    return Sigma

def CalculaDetJ(ne,XY,IJ):
    DetJ = np.zeros((ne))
    for e in range(ne):               
        # Recupera as coordenadas do elemento:
        X,Y = MontaXY(e,IJ,XY)
        
        dJ = np.zeros(4)
        for p in range(4):
            pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])
            
            r = pg[0,p]
            s = pg[1,p]
            
            #Calcula as derivadas da função de interpolação
            if X.size == 3:
                dNrs = MontadNTriangulo(r,s)
            elif X.size == 4:
                dNrs = MontadN(r,s)
            
            #Calcula a Matriz Jacobiana
            J = MontaJ(dNrs,X,Y)
            #Calcula determinante
            dJ[p] = np.linalg.det(J)

        DetJ[e] = min(dJ)
    
    return DetJ


def AplicaCCH(nn,na,AP,K_b,F_b):
    K = K_b.copy()
    F = F_b.copy()
    for l in range(na):
        #Nó
        no = AP[l,0]
        #gll
        gll = AP[l,1]
        
        #grau de liberdade global
        glg = 2*(no-1)+gll
        
        #Zera linha e coluna da Matriz
        for i in range(2*nn):
            K[i,glg-1] = 0
            K[glg-1,i] = 0
        #Coloca 1 na diagonal
        K[glg-1,glg-1] = 1
        
        F[glg-1] = 0
        
    return K,F

def AplicaCCLagrange(nn,nRL,RL,K_b,F_b):
    K0 = K_b.copy()
    K_aumentada = np.zeros((2*nn+nRL,2*nn+nRL))
    F0 = F_b.copy()
    F_aumentada = np.zeros(3*nn+nRL)

    K_aumentada[:2*nn,:2*nn] = K0
    F_aumentada[:2*nn] = F0
        
    for l in range(nRL):
        #Nó
        restricao=RL[l]
        
        if len(restricao) == 3:
            no = restricao[0]
            #gll
            gll = restricao[1]
            glg = 2*(no-1) + gll
            valor = restricao[2]
            K_aumentada[glg-1,2*nn+l] = 1
            K_aumentada[2*nn+l,glg-1] = 1
            F_aumentada[2*nn+l] = valor
            
        elif len(restricao) ==5:
            no1 = restricao[0]
            gll1 = restricao[1]
            no2 = restricao[2]
            gll2 = restricao[3]
            tipo = restricao[4]
            
            glg1 = 2*(no1-1)+gll1
            glg2 = 2*(no2-1)+gll2
            K_aumentada[2*nn+l,glg1-1] = 1
            K_aumentada[2*nn+l,glg2-1] = -tipo
            
            K_aumentada[glg1-1,2*nn+l] = 1
            K_aumentada[glg2-1,2*nn+l] = -tipo
            F_aumentada[2*nn+l] = 0
        
    return K_aumentada,F_aumentada


def TensaoVonMises2D(ne,sigma):
    TVM = np.zeros(ne)
    for i in range(ne):
        sig_xx = sigma[3*(i)]
        sig_yy = sigma[3*(i)+1]
        tau_xy = sigma[3*(i)+2]
        TVM[i] = np.sqrt(sig_xx**2 - sig_xx*sig_yy + sig_yy**2 + 3*tau_xy**2)
    
    return TVM


def RigidezGlobalBolha(nn,ne,MAT,ESP,XY,IJ):
    K = np.zeros((2*nn,2*nn))
    
    #Laço nos elementos
    for e in range(ne):
        Ee = MAT[e,0]
        ve = MAT[e,1]
        te = ESP[e]
        
        # Recupera as coordenadas do elemento:
        X,Y = MontaXY(e,IJ,XY)
        
        #Monta Ke
        Ke = MontaKeBolha(X,Y,Ee,ve,te)
        
        #Determina os gls globais do elemento:
        g = Montagls(e,IJ)
        
        for i in range(8):
            for j in range(8):
                K[int(g[i]),int(g[j])] += Ke[i,j]
                
        
    return K

def MontaKeBolha(X,Y,E,v,te):
    # Aloca Matriz
    Ke = np.zeros((8,8))
    Ke_ = np.zeros((12,12))
    Kuu = np.zeros((8,8))
    Kalfalfa = np.zeros((4,4))
    Kualfa = np.zeros((4,8))
    #Matriz Constitutiva do elemento
    
    Cv = MontaCEPT(E,v)
    
    #Define os pontos de Gauss
    
    pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])
    
    #Laço pelos ponos de Gauss
    for p in range(4):
        r = pg[0,p]
        s = pg[1,p]
        
        #Calcula as derivadas da função de interpolação
        dNrs = MontadNBolha(r,s)
        
        #Calcula a Matriz Jacobiana
        J = MontaJ(dNrs,X,Y)
        
        #Calcula determinante
        dJ = np.linalg.det(J)
        
        #Mapeia as derivadas para xy
        dNxy = CorrigedN(dNrs,J)
        
        #Monta B
        B = MontaBBolha(dNxy)

        #Acumula na matriz de rigidez do elemento

        Ke_ += ((np.transpose(B))@Cv@B)*te*dJ
    
    Kuu = Ke_[:8,:8]
    Kalfalfa = Ke_[8:,8:]
    Kualfa = Ke_[:8,8:]
    Kalfau = Kualfa.transpose()
    
    Ke = Kuu-(Kualfa@(np.linalg.inv(Kalfalfa))@Kalfau)
    
    
    return Ke


def MontadNBolha(r,s):
    #Matriz com as derivadas em relação a r (primeira linha) e s(segunda linha)
    
    dNrs = (0.25)*np.array([[s-1,1-s,1+s,-(1+s),4*(-2*r),0],[r-1,-(1+r),1+r,1-r,0,4*(-2*s)]])
    
    return dNrs


def MontaBBolha(dNxy):
    B = np.zeros((3,12))
    
    c = 0
    
    for i in range(6):
        B[0,c] += dNxy[0,i]
        B[1,c+1] += dNxy[1,i]
        B[2,c] += B[1,c+1]
        B[2,c+1] += B[0,c]
        
        c += 2
    return B


def CalculaTensaoMalhaBolha(nn,ne,MAT,ESP,XY,IJ,U):  
    
    Sigma = np.zeros((3*ne,1))
    
    for e in range(ne):
        Ue = np.zeros((8,1))
        Ee = MAT[e,0]
        ve = MAT[e,1]
        te = ESP[e]
        
        # Recupera as coordenadas do elemento:
        X,Y = MontaXY(e,IJ,XY)
        # Determina os gdls globais do elemento
        g = Montagls(e+1,IJ)
        
        for i in range(8):
            Ue[i] = U[int(g[i])]
            
        
        sigmae = CalculaTensaoElementoBolha(0,0,X,Y,Ee,ve,Ue,te)
        for j in range(3):
            Sigma[3*e+j] = sigmae[j]
        
    
    return Sigma

def CalculaTensaoElementoBolha(r,s,X,Y,E,v,Ue,te):
    #U(Ue+alfae)
    U = np.zeros((12,1))
    
    #Matriz constitutiva do elemento
    Cv = MontaCEPT(E,v)
    #Calcula as derivadas das funções de interpolação
    dNrs = MontadNBolha(r,s)
    #Calcula Jacobiana
    J = MontaJ(dNrs,X,Y)    
    #Mapeia as derivadas para xy
    dNxy = CorrigedN(dNrs,J)
    
    #Monta B
    B = MontaBBolha(dNxy)
    
    
    # Aloca Matriz
    Ke_ = np.zeros((12,12))
    Kalfalfa = np.zeros((4,4))
    Kualfa = np.zeros((4,8))
    #Matriz Constitutiva do elemento
    
    Cv = MontaCEPT(E,v)
    
    #Define os pontos de Gauss
    
    pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])
    
    #Laço pelos ponos de Gauss
    for p in range(4):
        r = pg[0,p]
        s = pg[1,p]
        
        #Calcula as derivadas da função de interpolação
        dNrs = MontadNBolha(r,s)
        
        #Calcula a Matriz Jacobiana
        J = MontaJ(dNrs,X,Y)
        
        #Calcula determinante
        dJ = np.linalg.det(J)
        
        #Mapeia as derivadas para xy
        dNxy = CorrigedN(dNrs,J)
        
        #Monta B
        B = MontaBBolha(dNxy)

        #Acumula na matriz de rigidez do elemento

        Ke_ += ((np.transpose(B))@Cv@B)*te*dJ
    
    Kalfalfa = Ke_[8:,8:]
    Kualfa = Ke_[8:,:8]
    
    alfae = -(np.linalg.inv(Kalfalfa))@Kualfa@Ue
    
    U[:8] = Ue
    for i in range(4):
        U[8+i] = alfae[i]
    
    sigmae = Cv@B@U
    
    return sigmae




def MontaKeTriangulo(X,Y,E,v,te):
    # Aloca Matriz
    Ke = np.zeros((6,6))
    
    #Matriz Constitutiva do elemento
    
    Cv = MontaCEPT(E,v) 
    
    #Define os pontos de Gauss
    
    pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])
    
    #Laço pelos ponos de Gauss
    for p in range(4):
        r = pg[0,p]
        s = pg[1,p]
        
        #Calcula as derivadas da função de interpolação
        dNrs = MontadNTriangulo(r,s)
        
        #Calcula a Matriz Jacobiana
        J = MontaJ(dNrs,X,Y)
        
        #Calcula determinante
        dJ = np.linalg.det(J)
        
        #Monta B
        B = MontaBTriagulo(X,Y)

        #Acumula na matriz de rigidez do elemento

        Ke += ((np.transpose(B))@Cv@B)*te*dJ
    return Ke
    
    
    
    
def MontadNTriangulo(r,s):
    #Matriz com as derivadas em relação a r (primeira linha) e s(segunda linha)
    
    dNrs = (0.25)*np.array([[s-1,1-s,0],[r-1,-(1+r),2]])
    
    return dNrs


def MontaBTriagulo(X,Y):

    x1 = X[0]
    x2 = X[1]
    x3 = X[2]
    y1 = Y[0]
    y2 = Y[1]
    y3 = Y[2]
    
    deltae = (x2-x1)*y3 + (x1-x3)*y2 + (x3-x2)*y1
    
    Bt = (1/deltae)*np.array([[y2-y3,0,y3-y1,0,y1-y2,0],[0,x3-x2,0,x1-x3,0,x2-x1],[x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]])
    
    return Bt


def FacesTriangulo(face):
    #Retorna dois índices para o cálculo dos vetores locais da face
    
    if face == 1 :
        p1 = 1
        p2 = 2
    elif face == 2:
        p1 =2
        p2 =3    
    elif face == 3:
        p1 =3
        p2=1
    return p1,p2

def FacesnvTriangulo(p1,p2,e,IJ,XY):
    #Recupera os nós da face
    n1 = IJ[e-1,p1-1]
    n2 = IJ[e-1,p2-1]
    
    #Recupera as coordenadas
    x1 = XY[n1-1,0]
    y1 = XY[n1-1,1]
    x2 = XY[n2-1,0]
    y2 = XY[n2-1,1]
    
    #Deltas
    deltax = x2-x1
    deltay = y2-y1
    
    L = (deltax**2 + deltay**2)**(0.5)
    #Determinante do Jacobiano
    dJ = L/2
    
    # Vetor tangente à face
    v = (1/L) * np.array([[deltax],[deltay]])
    # Vetor normal à face
    n = (1/L) * np.array([[deltay],[-deltax]])
    
    return dJ,n,v

def FtTriangulo(dJ,d,p1,p2,valor,te):
    #Pontos de gauss
    
    pg = (1/(3**(0.5))) * np.array([[-1],[1]])
    
    #Aloca o vetor
    
    F = np.zeros((6,1))
    
    #Laço nos pontos de Gauss:
    
    for i in range(2):
        ζ = pg[i]
        #Monta matriz N
        N =MontaNTriangulo(ζ,p1,p2)
        
        F = F+ valor*(N.transpose()@d)*dJ*te
        
    return F

def MontaNTriangulo(ζ,p1,p2):
    #Aloca N
    N = np.zeros((2,6))
    
    #Funções de interpolação não-nulas
    
    N1 = (1-ζ)/2
    N2 = (1+ζ)/2
    
    #Posiciona na Matriz N
    N[0,2*p1-2] = N1
    N[1,2*p1-1] = N1
    N[0,2*p2-2] = N2
    N[1,2*p2-1] = N2
    
    return N


def CalculaTensaoElementoTriangulo(r,s,X,Y,E,v,Ue):
    #Matriz constitutiva do elemento
    Cv = MontaCEPT(E,v)
    
    #Monta B
    B = MontaBTriagulo(X,Y)
    
    sigmae = Cv@B@Ue
    
    return sigmae    