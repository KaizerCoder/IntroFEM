import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def MontaLeθ(XY, ne, IJ):
    #Aloca Vetor de saída
    L = np.zeros(ne)
    θ = np.zeros(ne)
    
    for e in range(0, ne):
        #Nós do elemento
        n1 = IJ[e,0]
        n2 = IJ[e,1]
        
        #Coordenadas do nós
        x1 = XY[n1-1,0]
        y1 = XY[n1-1,1]
        
        x2 = XY[n2-1,0]
        y2 = XY[n2-1,1]
        l = ((x2-x1)**2 + (y2-y1)**2)**(0.5)
        L[e] = l
        θ[e] = np.arctan2(y2-y1,x2-x1)
        
    
    return L,θ
        
def MontaKGlobal(nn,ne,E,A,I,L,θ,IJ):
    K = np.zeros((3*nn,3*nn))
    
    for e in range(ne):
    
        E_e=E[e]
        A_e = A[e]
        L_e = L[e]
        I_e = I[e]
        θ_e = θ[e]
        
        #Nós do elemento
        n1 = IJ[e,0]
        n2 = IJ[e,1]
        
        #Rigidez local
        K_e = np.zeros((6,6))
        
        K_e = RigidezPP(E_e,A_e,I_e,L_e)
        
        #Matriz de Rotação do elemento
        T_e = MatrizTPP(θ_e)
        
        # Transforma a Rigidez do elemeto
        K_e_g = T_e.transpose()@K_e@T_e
        
        #Determina os gls do elemento
        
        g_e = glsPP(n1,n2)
        
        #Sobrepõe a matriz do elemento na rigidez global
        
        for i in range(6):
            for j in range(6):
                a = g_e[i]-1
                b = g_e[j]-1
                K[g_e[i]-1,g_e[j]-1] += K_e_g[i,j]
                
    return K

def RigidezPP(E,A,I,L):
    
    K_l = np.array([
    [ E*A/L,          0,             0,          -E*A/L,            0,              0],
    [ 0,      12*E*I/L**3,    6*E*I/L**2,            0,   -12*E*I/L**3,    6*E*I/L**2],
    [ 0,       6*E*I/L**2,     4*E*I/L,             0,    -6*E*I/L**2,     2*E*I/L],
    [-E*A/L,          0,             0,           E*A/L,             0,              0],
    [ 0,     -12*E*I/L**3,   -6*E*I/L**2,           0,    12*E*I/L**3,     -6*E*I/L**2],
    [ 0,      6*E*I/L**2,    2*E*I/L,              0,     -6*E*I/L**2,     4*E*I/L]
])   
    # print(check_symmetric(K_l))

    
    return K_l

def MatrizTPP(θ):
    s = np.sin(θ)
    c = np.cos(θ)
    
    T = np.array([
    [ c,  s, 0,  0,  0, 0],
    [-s,  c, 0,  0,  0, 0],
    [ 0,  0, 1,  0,  0, 0],
    [ 0,  0, 0,  c,  s, 0],
    [ 0,  0, 0, -s,  c, 0],
    [ 0,  0, 0,  0,  0, 1]
])
    return T

def glsPP(n1,n2):
    
    g = np.array([  3*(n1-1)+1,
                    3*(n1-1)+2,
                    3*(n1-1)+3,
                    3*(n2-1)+1,
                    3*(n2-1)+2,
                    3*(n2-1)+3])
    
    return g

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
    

def MontaFGlobal(nn,nf,FC,L,θ,IJ,nq,FQ):
    F = np.zeros(3*nn)
    
    Fc = FCGlobal(nn,nf,FC)
    Fq = FQGlobal(nn,L,θ,IJ,nq,FQ)
    
    F = Fc+Fq
    
    return F


def FCGlobal(nn,nf,FC):
    Fc = np.zeros(3*nn)
    
    for i in range(nf):
        #Nó
        no = FC[i,0]
        #Grau de liberdade
        gll = FC[i,1]
        #Valor da Força
        valor = FC[i,2]
        #Grau de liberdade global
        glg = 3*(no-1)+gll
        
        Fc[glg-1] += valor
    return Fc

def FQGlobal(nn,L,θ,IJ,nq,FQ):
    Fq = np.zeros(3*nn)
    
    for i in range(nq):
        #elemento
        e = FQ[i,0]
        #Valor da Força
        q1 = FQ[i,1]
        q2 = FQ[i,2]
        
        #Ângulo e Comprimeno
        Le = L[e-1]
        θe = θ[e-1]
        
        Qe = FQPP(Le,q1,q2)
        Te = MatrizTPP(θe)
        
        Qeg = Te.transpose()@VetorColuna(Qe) 
        
        n1 = IJ[e-1,0]
        n2 = IJ[e-1,1]
        
        g= glsPP(n1,n2)
        
        for j in range(6):
            Fq[g[j]-1] += Qeg[j] 
        
    return Fq


def FQPP(L,q1,q2):
    
    Qel = np.array([0,
                    (3*L*q2+7*L*q1)/20,
                    (2*(L**2)*q2+3*(L**2)*q1)/60,
                    0,
                    (7*L*q2+3*L*q1)/20,
                    -(3*(L**2)*q2+2*(L**2)*q1)/60,
                    ])
    
    return Qel

def EsforcosInternos(ne, E, A, I, L, θ, IJ, U):
    esforços = []

    for e in range(ne):
        # Parâmetros do elemento
        E_e = E[e]
        A_e = A[e]
        I_e = I[e]
        L_e = L[e]
        θ_e = θ[e]

        # Rigidez local e matriz de rotação
        K_e = RigidezPP(E_e, A_e, I_e, L_e)
        T_e = MatrizTPP(θ_e)

        # GLS do elemento e deslocamentos globais correspondentes
        n1 = IJ[e, 0]
        n2 = IJ[e, 1]
        g = glsPP(n1, n2)
        Ue_global = U[g - 1] 

        # Transforma para vetor coluna local
        Ue_global = np.array(Ue_global).reshape((6, 1))
        Ue_local = T_e @ Ue_global

        # Calcula os esforços internos no sistema local
        Fe_local = -K_e @ Ue_local
                
        esforços.append(Fe_local.flatten())  # Pode retornar como lista plana

    return esforços

def AplicaCCLagrange(nn,nRL,RL,K_b,F_b):
    K0 = K_b.copy()
    K_aumentada = np.zeros((3*nn+nRL,3*nn+nRL))
    F0 = F_b.copy()
    F_aumentada = np.zeros(3*nn+nRL)

    K_aumentada[:3*nn,:3*nn] = K0
    F_aumentada[:3*nn] = F0
        
    for l in range(nRL):
        #Nó
        restricao=RL[l]
        
        if len(restricao) == 3:
            no = restricao[0]
            #gll
            gll = restricao[1]
            glg = 3*(no-1) + gll
            valor = restricao[2]
            K_aumentada[glg-1,3*nn+l] = 1
            K_aumentada[3*nn+l,glg-1] = 1
            F_aumentada[3*nn+l] = valor
            
        elif len(restricao) ==5:
            no1 = restricao[0]
            gll1 = restricao[1]
            no2 = restricao[2]
            gll2 = restricao[3]
            tipo = restricao[4]
            
            glg1 = 3*(no1-1)+gll1
            glg2 = 3*(no2-1)+gll2
            K_aumentada[3*nn+l,glg1-1] = 1
            K_aumentada[3*nn+l,glg2-1] = -tipo
            
            K_aumentada[glg1-1,3*nn+l] = 1
            K_aumentada[glg2-1,3*nn+l] = -tipo
            F_aumentada[3*nn+l] = 0
        
    return K_aumentada,F_aumentada


def VetorColuna(vec):
    return np.reshape(vec, (len(vec), 1))

def check_symmetric(a, tol=1e-8):
    return np.all(np.abs(a-a.T) < tol)
    