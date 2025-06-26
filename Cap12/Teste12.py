import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import PorticoPlano as PP


#Entrada de Dados 

#Número de nós
nn = 4

#Coordenadas
# XY = np.array('0 0;0 1; 2 1')
XY = np.array([[0,0],[0,1],[2,1],[2,0.5]])

#Número de elementos
ne = 3

#Conectividade
IJ = np.array([[1,2],[2,3],[3,4]])

#Parâmetros dos elementos
E = 210e9*np.array([1,1,1])

A = 0.01*np.array([1,1,1])

I = 0.0002*np.array([1,1,1])

#Forças concentradas
nf = 4

FC  = np.array([[2,1,3000],[2,2,-3000],[2,3,100],[3,1,5000]])

nq = 0
FQ = np.zeros(0)

# nRL = 5
#Exemplos de CCs:
#e1 gdl 1 = 0 -> [1,1,0]
#e3 gdl 2 = 0.005 -> [3,2,0.005]
# e2 gdl2 = e3 gdl2 -> [2,2,3,2,1]
# e2 gdl2 = -e3 gdl2(rótula) -> [2,2,3,2,-1]
# e2 gdl2 = 3*e3 gdl2(rótula) -> [2,2,3,2,3]


# RL = np.array([[1,1,0],[1,2,0],[4,1,0],[4,2,0],[2,2,3,2]],dtype=object)
nRL =4
RL = np.array([[1,1,0],[1,2,0],[4,1,0],[4,2,0]])

#Pré Processamento 
L,θ = PP.MontaLeθ(XY, ne, IJ)

K_b = PP.MontaKGlobal(nn,ne,E,A,I,L,θ,IJ)

F = PP.MontaFGlobal(nn,nf,FC,L,θ,IJ,nq,FQ)

K_modL,F_modL = PP.AplicaCCLagrange(nn,nRL,RL,K_b,F) 

U_L = np.linalg.solve(K_modL,F_modL)

U = U_L[:3*nn]

sp.pprint(U)

esforços = PP.EsforcosInternos(ne, E, A, I, L, θ, IJ, U)
pass

