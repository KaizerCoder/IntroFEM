import numpy as np

#Rigidez Local de Elemento:

#Montagem da matriz constitutiva EPT de um elemento com E e v para caso 3D
def MontaCEPT3D(E,v):
    c = E/((1+v)*(1-2*v))
    
    Cv = c*np.array([
        [1-v,v,v,0,0,0],
        [v,1-v,v,0,0,0],
        [v,v,1-v,0,0,0],
        [0,0,0,(1-2*v)/2,0,0],
        [0,0,0,0,(1-2*v)/2,0],
        [0,0,0,0,0,(1-2*v)/2],
        ])
    
    return Cv

def MontadN(r,s,t):
    #Matriz com as derivadas em relação a r (primeira linha) e s(segunda linha) e t (terceira linha)
    
    dNdr = (1/8) * np.array([
        -(1 - s)*(1 - t),
         (1 - s)*(1 - t),
         (1 + s)*(1 - t),
        -(1 + s)*(1 - t),
        -(1 - s)*(1 + t),
         (1 - s)*(1 + t),
         (1 + s)*(1 + t),
        -(1 + s)*(1 + t)
    ])

    dNds = (1/8) * np.array([
        -(1 - r)*(1 - t),
        -(1 + r)*(1 - t),
         (1 + r)*(1 - t),
         (1 - r)*(1 - t),
        -(1 - r)*(1 + t),
        -(1 + r)*(1 + t),
         (1 + r)*(1 + t),
         (1 - r)*(1 + t)
    ])

    dNdt = (1/8) * np.array([
        -(1 - r)*(1 - s),
        -(1 + r)*(1 - s),
        -(1 + r)*(1 + s),
        -(1 - r)*(1 + s),
         (1 - r)*(1 - s),
         (1 + r)*(1 - s),
         (1 + r)*(1 + s),
         (1 - r)*(1 + s)
    ])

    dNrst = np.vstack((dNdr, dNds, dNdt))
    
    return dNrst

def MontaJ(dNrst,X,Y,Z):
    
    J = np.zeros((3,3))
    
    for i in range(X.size):
        J[0, 0] += dNrst[0, i] * X[i]
        J[0, 1] += dNrst[0, i] * Y[i]
        J[0, 2] += dNrst[0, i] * Z[i]
        
        J[1, 0] += dNrst[1, i] * X[i]
        J[1, 1] += dNrst[1, i] * Y[i]
        J[1, 2] += dNrst[1, i] * Z[i]
        
        J[2, 0] += dNrst[2, i] * X[i]
        J[2, 1] += dNrst[2, i] * Y[i]
        J[2, 2] += dNrst[2, i] * Z[i]
    
    return J

def CorrigedN(dNrs,J):

    dNxyz = (np.linalg.inv(J))@dNrs
    
    return dNxyz

def MontaB(dNxyz):
    B = np.zeros((6,24))
        
    for i in range(8):
        idx = 3 * i
        dN_dx = dNxyz[0, i]
        dN_dy = dNxyz[1, i]
        dN_dz = dNxyz[2, i]
        
        B[0, idx    ] = dN_dx
        B[1, idx + 1] = dN_dy
        B[2, idx + 2] = dN_dz
        B[3, idx    ] = dN_dy
        B[3, idx + 1] = dN_dx
        B[4, idx + 1] = dN_dz
        B[4, idx + 2] = dN_dy
        B[5, idx    ] = dN_dz
        B[5, idx + 2] = dN_dx
        
    return B

#Monta Ke

def MontaKe(X,Y,Z,E,v):
    # Aloca Matriz
    Ke = np.zeros((24,24))
    #Matriz Constitutiva do elemento
    
    Cv = MontaCEPT3D(E,v) #
    
    #Define os pontos de Gauss
    
    g = 1.0 / np.sqrt(3)
    pontos_gauss = [(-g, -g, -g), (g, -g, -g), (g, g, -g), (-g, g, -g),
                    (-g, -g, g),  (g, -g, g),  (g, g, g),  (-g, g, g)]    
    #Laço pelos ponos de Gauss
    for (r,s,t) in pontos_gauss:
        
        #Calcula as derivadas da função de interpolação
        dNrst = MontadN(r,s,t)
        
        #Calcula a Matriz Jacobiana
        J = MontaJ(dNrst,X,Y,Z)
        
        #Calcula determinante
        dJ = np.linalg.det(J)
        
        #Mapeia as derivadas para xy
        dNxyz = CorrigedN(dNrst,J)
        
        #Monta B
        B = MontaB(dNxyz)

        #Acumula na matriz de rigidez do elemento

        Ke += ((np.transpose(B))@Cv@B)*dJ
    return Ke

def MontaXYZ(e,IJ,XYZ):
    nos = IJ[e-1]
    
    X = np.zeros((nos.size))
    Y = np.zeros((nos.size))
    Z = np.zeros((nos.size))
    
    #Preenche os vetores
    
    for i in range(nos.size):
        #Nó
        no = IJ[e-1,i]
        #Coordenadas
        X[i] = XYZ[no-1,0]
        Y[i] = XYZ[no-1,1]
        Z[i] = XYZ[no-1,2]
        
    return X,Y,Z

def Montagls(e,IJ):
    #Inicializa vetor
    nos = IJ[e-1]

    gls = np.zeros((3*nos.size,1))      
    c = 0
    
    for i in range(nos.size):
        #Nó
        no = nos[i]
        #Laço pelos graus de liberdade dos nós:
        for j in range(3):
            gls[c] =  int(3*(no-1)+j)
            c+=1
    
    return gls
        
    
def RigidezGlobal(nn,ne,MAT,XYZ,IJ,bolha = False):
    K = np.zeros((3*nn,3*nn))
    
    #Laço nos elementos
    for e in range(ne):
        Ee = MAT[e,0]
        ve = MAT[e,1]
        
        #Determina os gls globais do elemento:
        g = Montagls(e,IJ)
        
        # Recupera as coordenadas do elemento:
        X,Y,Z = MontaXYZ(e,IJ,XYZ)
        
        #Monta Ke
        if g.size == 24:
            if bolha:
                Ke = MontaKeBolha(X,Y,Z,Ee,ve)
            else:
                Ke = MontaKe(X,Y,Z,Ee,ve)
        # elif g.size == 6:
        #     Ke = MontaKeTriangulo(X,Y,Ee,ve)
        
        for i in range(g.size):
                for j in range(g.size):
                    K[int(g[i]),int(g[j])] += Ke[i,j]
    return K

def Faces(face):
    #Retorna as propriedades da Quatro índices para o cálculo dos vetores locais da face e o indice a
    
    #1234 e t=-1 fixo
    if face == 1 :
        p1 = 1
        p2 = 2
        p3 = 3
        p4 = 4
        a = 1
    #5678 e t=1 fixo
    elif face == 2:
        p1 = 5
        p2 = 6   
        p3 = 7
        p4 = 8
        a=0
    #1265 e r=-1 fixo
    elif face == 3:
        p1 = 1
        p2 = 2
        p3 = 6
        p4 = 5
        a=0
    #2376 e r=1 fixo
    elif face == 4:
        p1 = 2
        p2 = 3
        p3 = 7
        p4 = 6
        a=0
    #4378 e s=-1 fixo
    elif face == 5:
        p1 = 4
        p2 = 3
        p3 = 7
        p4 = 8
        a=1
    #1487 e s=1 fixo
    elif face == 6:
        p1 = 1
        p2 = 4
        p3 = 8
        p4 = 7
        a=0
    return p1,p2,p3,p4,a

def MontaN(ζ,η,p1,p2,p3,p4):
    #Aloca N
    N = np.zeros((3,24))
    
    #Funções de interpolação não-nulas
    
    N1 = 0.25*(1-ζ)*(1-η)
    N2 = 0.25*(1+ζ)*(1-η)
    N3 = 0.25*(1+ζ)*(1+η)
    N4 = 0.25*(1-ζ)*(1+η)
    
    #Posiciona na Matriz N
    for i, Ns in zip([p1, p2, p3, p4], [N1, N2, N3, N4]):
        c = 3 * (i - 1)
        N[0,c] = Ns
        N[1,c+ 1] = Ns
        N[2,c+2] = Ns
    
    return N

def Facesnv(p1, p2, p4,a, e, IJ, XYZ):
    #Recupera os nós da face
    n1 = IJ[e-1,p1-1]
    n2 = IJ[e-1,p2-1]
    n4 = IJ[e-1,p4-1]

    
    #Recupera as coordenadas
    P1 = XYZ[n1-1]
    P2 = XYZ[n2-1]
    P4 = XYZ[n4-1]
    
    #Vetores tangentes da face
    vζ = (P2 - P1)/2
    vη = (P4 - P1)/2
    
    n = ((-1)**a)*np.cross(vζ,vη)
    
    
    dA = np.linalg.norm(n)
    n = n/np.linalg.norm(n)
    vζ = vζ/np.linalg.norm(vζ)
    vη = vη/np.linalg.norm(vη)

           
    return n,vζ,vη,dA

def FtFaceH8(dA,d,p1,p2,p3,p4,val):   
    #Aloca o vetor
    
    F = np.zeros((24))
    
    #Laço nos pontos de Gauss:
    gauss = [-1 / np.sqrt(3), 1 / np.sqrt(3)]
    for ζ in gauss:
        for η in gauss:
            N = MontaN(ζ,η,p1,p2,p3,p4)
            F = F+ val*(N.transpose()@d)*dA
        
    return F

#Montagem do vetor global de forças distribuídas nas faces dos elementos e diretamente nos nós
def ForcapGlobal(nn,XYZ,IJ,nc,P):
    F = np.zeros((3*nn))
    
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
            F[3*(no_global-1)+(gdl_local-1)] += val    
        
        #Força numa face
        else:
            #Determina os gls do elemento
            g =  Montagls(ele,IJ)
            
            
            if g.size == 24:
                #Propriedades da face
                p1,p2,p3,p4,a = Faces(face)
                
                #Vetores da face
                n,vζ,vη,dA = Facesnv(p1,p2,p4,a,ele,IJ,XYZ)
                
                if dire == 1:
                    d = n
                elif dire == 2:
                    d = vζ
                elif dire == 3:
                    d = vη
                    
                #Vetor de Força
                Fe = FtFaceH8(dA,d,p1,p2,p3,p4,val)
                
            # elif g.size == 6:
            #     p1,p2 = FacesTriangulo(face)
            #     #Vetores da face
            #     dJ,n,v = FacesnvTriangulo(p1,p2,ele,IJ,XYZ)
                
            #     if dire == 1:
            #         d = n
            #     elif dire == 2:
            #         d = v
            #     Fe = FtTriangulo(dJ,d,p1,p2,val,te)
            
            
            for j in range(g.size):
                    F[int(g[j])] += Fe[j]
    return F
        
#Retorna o vetor com as tensões em um ponto (r,s) de um elemento bilinear isoparamétrico de 4 nós.

def CalculaTensaoElemento(r,s,t,X,Y,Z,E,v,Ue):
    #Matriz constitutiva do elemento
    Cv = MontaCEPT3D(E,v)
    #Calcula as derivadas das funções de interpolação
    dNrst = MontadN(r,s,t)
    #Calcula Jacobiana
    J = MontaJ(dNrst,X,Y,Z)    
    #Mapeia as derivadas para xy
    dNxyz = CorrigedN(dNrst,J)
    
    #Monta B
    B = MontaB(dNxyz)
    
    sigmae = Cv@B@Ue
    
    return sigmae
    

def CalculaTensaoMalha(nn,ne,MAT,XYZ,IJ,U,bolha=False):  
    
    Sigma = np.zeros(6*ne)
    
    for e in range(ne):
        
        Ee = MAT[e,0]
        ve = MAT[e,1]
                
        # Recupera as coordenadas do elemento:
        X,Y,Z = MontaXYZ(e,IJ,XYZ)
        # Determina os gdls globais do elemento
        g = Montagls(e+1,IJ)
        
        Ue = np.zeros((g.size,1))

        for i in range(g.size):
            Ue[i] = U[int(g[i])]
        
        if g.size == 24:
            if bolha:
                sigmae = CalculaTensaoElementoBolha(0,0,0,X,Y,Z,Ee,ve,Ue)
            else:
                sigmae = CalculaTensaoElemento(0,0,0,X,Y,Z,Ee,ve,Ue)
        # elif g.size == 6:
        #     sigmae = CalculaTensaoElementoTriangulo(0,0,X,Y,Ee,ve,Ue)
        for j in range(6):
            Sigma[6*e+j] = sigmae[j]
        
    
    return Sigma

# def CalculaDetJ(ne,XY,IJ):
#     DetJ = np.zeros((ne))
#     for e in range(ne):               
#         # Recupera as coordenadas do elemento:
#         X,Y = MontaXY(e,IJ,XY)
        
#         dJ = np.zeros(4)
#         for p in range(4):
#             pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])
            
#             r = pg[0,p]
#             s = pg[1,p]
            
#             #Calcula as derivadas da função de interpolação
#             if X.size == 3:
#                 dNrs = MontadNTriangulo(r,s)
#             elif X.size == 4:
#                 dNrs = MontadN(r,s)
            
#             #Calcula a Matriz Jacobiana
#             J = MontaJ(dNrs,X,Y)
#             #Calcula determinante
#             dJ[p] = np.linalg.det(J)

#         DetJ[e] = min(dJ)
    
#     return DetJ


def AplicaCCH(nn,na,AP,K_b,F_b):
    K = K_b.copy()
    F = F_b.copy()
    for l in range(na):
        #Nó
        no = AP[l,0]
        #gll
        gll = AP[l,1]
        
        #grau de liberdade global
        glg = 3*(no-1)+gll
        
        #Zera linha e coluna da Matriz
        for i in range(3*nn):
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


def TensaoVonMises3D(ne,sigma):
    TVM = np.zeros(ne)
    for i in range(ne):
        sig_xx = sigma[6 * i]
        sig_yy = sigma[6 * i + 1]
        sig_zz = sigma[6 * i + 2]
        tau_xy = sigma[6 * i + 3]
        tau_yz = sigma[6 * i + 4]
        tau_zx = sigma[6 * i + 5]

        TVM[i] = np.sqrt(
            0.5 * (
                (sig_xx - sig_yy) ** 2 +
                (sig_yy - sig_zz) ** 2 +
                (sig_zz - sig_xx) ** 2
            ) + 3 * (
                tau_xy ** 2 + tau_yz ** 2 + tau_zx ** 2
            )
        )
    
    return TVM

def MontaKeBolha(X,Y,Z,E,v):
    # Aloca Matriz
    Ke = np.zeros((24,24))
    Ke_ = np.zeros((33,33))
    Kuu = np.zeros((24,24))
    Kalfalfa = np.zeros((9,9))
    Kualfa = np.zeros((9,24))
    #Matriz Constitutiva do elemento
    
    Cv = MontaCEPT3D(E,v)
    
    #Laço pelos ponos de Gauss
    g = 1.0 / np.sqrt(3)
    pontos_gauss = [(-g, -g, -g), (g, -g, -g), (g, g, -g), (-g, g, -g),
                    (-g, -g, g),  (g, -g, g),  (g, g, g),  (-g, g, g)]    
    #Laço pelos pontos de Gauss
    for (r,s,t) in pontos_gauss:
        
        #Calcula as derivadas da função de interpolação
        dNrst = MontadNBolha(r,s,t)
        
        #Calcula a Matriz Jacobiana
        J = MontaJ(dNrst,X,Y,Z)
        
        #Calcula determinante
        dJ = np.linalg.det(J)
        
        #Mapeia as derivadas para xy
        dNxyz = CorrigedN(dNrst,J)
        
        #Monta B
        B = MontaBBolha(dNxyz)

        #Acumula na matriz de rigidez do elemento

        Ke_ += ((np.transpose(B))@Cv@B)*dJ
    
    Kuu = Ke_[:24,:24]
    Kalfalfa = Ke_[24:,24:]
    Kualfa = Ke_[:24,24:]
    Kalfau = Kualfa.transpose()
    
    Ke = Kuu-(Kualfa@(np.linalg.inv(Kalfalfa))@Kalfau)
    
    
    return Ke


def MontadNBolha(r,s,t):
    #Matriz com as derivadas em relação a r (primeira linha) e s(segunda linha)
    
    dNdr = (1/8) * np.array([
        -(1 - s)*(1 - t),
         (1 - s)*(1 - t),
         (1 + s)*(1 - t),
        -(1 + s)*(1 - t),
        -(1 - s)*(1 + t),
         (1 - s)*(1 + t),
         (1 + s)*(1 + t),
        -(1 + s)*(1 + t),
          8*(-2*r),
          0,
          0
    ])

    dNds = (1/8) * np.array([
        -(1 - r)*(1 - t),
        -(1 + r)*(1 - t),
         (1 + r)*(1 - t),
         (1 - r)*(1 - t),
        -(1 - r)*(1 + t),
        -(1 + r)*(1 + t),
         (1 + r)*(1 + t),
         (1 - r)*(1 + t),
         0,
         8*(-2*s),
         0
    ])

    dNdt = (1/8) * np.array([
        -(1 - r)*(1 - s),
        -(1 + r)*(1 - s),
        -(1 + r)*(1 + s),
        -(1 - r)*(1 + s),
         (1 - r)*(1 - s),
         (1 + r)*(1 - s),
         (1 + r)*(1 + s),
         (1 - r)*(1 + s),
         0,
         0,
         8*(-2*t)
    ])

    dNrst = np.vstack((dNdr, dNds, dNdt))
    
    return dNrst


def MontaBBolha(dNxyz):
    B = np.zeros((6,33))
        
    for i in range(11):
        idx = 3 * i
        dN_dx = dNxyz[0, i]
        dN_dy = dNxyz[1, i]
        dN_dz = dNxyz[2, i]
        
        B[0, idx    ] = dN_dx
        B[1, idx + 1] = dN_dy
        B[2, idx + 2] = dN_dz
        B[3, idx    ] = dN_dy
        B[3, idx + 1] = dN_dx
        B[4, idx + 1] = dN_dz
        B[4, idx + 2] = dN_dy
        B[5, idx    ] = dN_dz
        B[5, idx + 2] = dN_dx
    return B

def CalculaTensaoElementoBolha(r,s,t,X,Y,Z,E,v,Ue):
    #U(Ue+alfae)
    U = np.zeros((33,1))
    
    #Matriz constitutiva do elemento
    Cv = MontaCEPT3D(E,v)
    #Calcula as derivadas das funções de interpolação
    dNrst = MontadNBolha(r,s,t)
    #Calcula Jacobiana
    J = MontaJ(dNrst,X,Y,Z)    
    #Mapeia as derivadas para xy
    dNxyz = CorrigedN(dNrst,J)
    
    #Monta B
    B = MontaBBolha(dNxyz)
    
    
    # Aloca Matriz
    Ke_ = np.zeros((33,33))
    Kalfalfa = np.zeros((9,9))
    Kualfa = np.zeros((9,24))
    #Matriz Constitutiva do elemento
    
    Cv = MontaCEPT3D(E,v)
    
    #Define os pontos de Gauss
    g = 1.0 / np.sqrt(3)
    pontos_gauss = [(-g, -g, -g), (g, -g, -g), (g, g, -g), (-g, g, -g),
                    (-g, -g, g),  (g, -g, g),  (g, g, g),  (-g, g, g)]  
    
    #Laço pelos pontos de Gauss
    for (r,s,t) in pontos_gauss:
        
        #Calcula as derivadas da função de interpolação
        dNrst = MontadNBolha(r,s,t)
        
        #Calcula a Matriz Jacobiana
        J = MontaJ(dNrst,X,Y,Z)
        
        #Calcula determinante
        dJ = np.linalg.det(J)
        
        #Mapeia as derivadas para xy
        dNxyz = CorrigedN(dNrst,J)
        
        #Monta B
        B = MontaBBolha(dNxyz)

        #Acumula na matriz de rigidez do elemento

        Ke_ += ((np.transpose(B))@Cv@B)*dJ
    
    Kalfalfa = Ke_[24:,24:]
    Kualfa = Ke_[24:,:24]
    
    alfae = -(np.linalg.inv(Kalfalfa))@Kualfa@Ue
    
    U[:24] = Ue
    for i in range(4):
        U[24+i] = alfae[i]
    
    sigmae = Cv@B@U
    
    return sigmae




# def MontaKeTriangulo(X,Y,E,v,te):
#     # Aloca Matriz
#     Ke = np.zeros((6,6))
    
#     #Matriz Constitutiva do elemento
    
#     Cv = MontaCEPT(E,v) 
    
#     #Define os pontos de Gauss
    
#     pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])
    
#     #Laço pelos ponos de Gauss
#     for p in range(4):
#         r = pg[0,p]
#         s = pg[1,p]
        
#         #Calcula as derivadas da função de interpolação
#         dNrs = MontadNTriangulo(r,s)
        
#         #Calcula a Matriz Jacobiana
#         J = MontaJ(dNrs,X,Y)
        
#         #Calcula determinante
#         dJ = np.linalg.det(J)
        
#         #Monta B
#         B = MontaBTriagulo(X,Y)

#         #Acumula na matriz de rigidez do elemento

#         Ke += ((np.transpose(B))@Cv@B)*te*dJ
#     return Ke
    
    
    
    
# def MontadNTriangulo(r,s):
#     #Matriz com as derivadas em relação a r (primeira linha) e s(segunda linha)
    
#     dNrs = (0.25)*np.array([[s-1,1-s,0],[r-1,-(1+r),2]])
    
#     return dNrs


# def MontaBTriagulo(X,Y):

#     x1 = X[0]
#     x2 = X[1]
#     x3 = X[2]
#     y1 = Y[0]
#     y2 = Y[1]
#     y3 = Y[2]
    
#     deltae = (x2-x1)*y3 + (x1-x3)*y2 + (x3-x2)*y1
    
#     Bt = (1/deltae)*np.array([[y2-y3,0,y3-y1,0,y1-y2,0],[0,x3-x2,0,x1-x3,0,x2-x1],[x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]])
    
#     return Bt


# def FacesTriangulo(face):
#     #Retorna dois índices para o cálculo dos vetores locais da face
    
#     if face == 1 :
#         p1 = 1
#         p2 = 2
#     elif face == 2:
#         p1 =2
#         p2 =3    
#     elif face == 3:
#         p1 =3
#         p2=1
#     return p1,p2

# def FacesnvTriangulo(p1,p2,e,IJ,XY):
#     #Recupera os nós da face
#     n1 = IJ[e-1,p1-1]
#     n2 = IJ[e-1,p2-1]
    
#     #Recupera as coordenadas
#     x1 = XY[n1-1,0]
#     y1 = XY[n1-1,1]
#     x2 = XY[n2-1,0]
#     y2 = XY[n2-1,1]
    
#     #Deltas
#     deltax = x2-x1
#     deltay = y2-y1
    
#     L = (deltax**2 + deltay**2)**(0.5)
#     #Determinante do Jacobiano
#     dJ = L/2
    
#     # Vetor tangente à face
#     v = (1/L) * np.array([[deltax],[deltay]])
#     # Vetor normal à face
#     n = (1/L) * np.array([[deltay],[-deltax]])
    
#     return dJ,n,v

# def FtTriangulo(dJ,d,p1,p2,valor,te):
#     #Pontos de gauss
    
#     pg = (1/(3**(0.5))) * np.array([[-1],[1]])
    
#     #Aloca o vetor
    
#     F = np.zeros((6,1))
    
#     #Laço nos pontos de Gauss:
    
#     for i in range(2):
#         ζ = pg[i]
#         #Monta matriz N
#         N =MontaNTriangulo(ζ,p1,p2)
        
#         F = F+ valor*(N.transpose()@d)*dJ*te
        
#     return F

# def MontaNTriangulo(ζ,p1,p2):
#     #Aloca N
#     N = np.zeros((2,6))
    
#     #Funções de interpolação não-nulas
    
#     N1 = (1-ζ)/2
#     N2 = (1+ζ)/2
    
#     #Posiciona na Matriz N
#     N[0,2*p1-2] = N1
#     N[1,2*p1-1] = N1
#     N[0,2*p2-2] = N2
#     N[1,2*p2-1] = N2
    
#     return N


# def CalculaTensaoElementoTriangulo(r,s,X,Y,E,v,Ue):
#     #Matriz constitutiva do elemento
#     Cv = MontaCEPT(E,v)
    
#     #Monta B
#     B = MontaBTriagulo(X,Y)
    
#     sigmae = Cv@B@Ue
    
#     return sigmae    