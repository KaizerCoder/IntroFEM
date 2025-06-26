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
        
def MontaKGlobal(nn,ne,E,A,L,θ,IJ):
    K = np.zeros((2*nn,2*nn))
    
    for e in range(ne):
    
        E_e=E[e]
        A_e = A[e]
        L_e = L[e]
        θ_e = θ[e]
        
        #Rigidez local
        K_e = RigidezBarra(E_e,A_e,L_e)
        
        #Matriz de Rotação do elemento
        T_e = MatrizTBarra(θ_e)
        
        # Transforma a Rigidez do elemeto
        K_e_g = T_e.transpose()*K_e*T_e
        
        #Nós do elemento
        n1 = IJ[e,0]
        n2 = IJ[e,1]
        
        #Determina os gls do elemento
        
        g_e = glsBarra(n1,n2)
        pass
        
        #Sobrepõe a matriz do elemento na rigidez global
        
        for i in range(4):
            for j in range(4):
                a = g_e[i]-1
                b = g_e[j]-1
                K[g_e[i]-1,g_e[j]-1] += K_e_g[i,j]
                
    return K

def RigidezBarra(E,A,L):
    
    K_l = (E*A/L)*np.matrix('1 0 -1 0;0 0 0 0;-1 0 1 0;0 0 0 0')
    
    return K_l

def MatrizTBarra(θ):
    s = np.sin(θ)
    c = np.cos(θ)
    
    T = np.matrix([[c,s,0,0],[-s,c,0,0],[0,0,c,s],[0,0,-s,c]])    
    return T

def glsBarra(n1,n2):
    
    g = np.array([2*(n1-1)+1,2*(n1-1)+2,2*(n2-1)+1,2*(n2-1)+2])
    
    return g

def AplicaCCH(nn,na,AP,K,F):
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
    


def MontaFGlobal(nn,nf,FC):
    F = np.zeros(2*nn)
    
    for i in range(nf):
        #Nó
        no = FC[i,0]
        #Grau de liberdade
        gll = FC[i,1]
        #Valor da Força
        valor = FC[i,2]
        #Grau de liberdade global
        glg = 2*(no-1)+gll
        
        F[glg-1] += valor
    return F
        

def EsforcosInternos(ne,E,A,L,θ,IJ,U):
    
    N = np.zeros(ne)
    
    for e in range(ne):
        E_e=E[e]
        A_e = A[e]
        L_e = L[e]
        θ_e = θ[e]
        
        #Rigidez local
        K_e = RigidezBarra(E_e,A_e,L_e)
        
        #Matriz de Rotação do elemento
        T_e = MatrizTBarra(θ_e)
                
        #Nós do elemento
        n1 = IJ[e,0]
        n2 = IJ[e,1]
        
        #Determina os gls do elemento
        
        g_e = glsBarra(n1,n2)
        
        #Recupera os deslocamentos do elemento(global)
        U_e_g = (U[g_e-1])

        #Transforma em um vetor coluna
        U_e_g=VetorColuna(U_e_g)        
      
        #Transforma para sistema local
        U_e_l = T_e*U_e_g
        
        #Calcula o vetor de Forças internas do elemento
        F_e_l = K_e*U_e_l
        
        #Armazena o esforço interno
        N[e] = -F_e_l[0]
    return N
    

def VetorColuna(vec):
    return sp.Matrix(vec).T.T