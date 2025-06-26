import numpy as np
import BiLinearIso4Nos as BL4


#Entrada de Dados 

#Número de nós
nn = 6
XY = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.0],[0.0,0.1],[0.5,0.1],[1.0,0.1]])
ne = 2
IJ = np.array([[1,2,5,4],[2,3,6,5]])

MAT = np.array([[1e9,0.3],[1e9,0.3],[1e9,0.3]])

ESP = np.array([[0.01],[0.01],[0.01]])

na = 5
AP = np.array([[1,1,0],[1,2,0],[2,2,0],[3,2,0],[4,1,0]])

nc=1
P = np.array([[2,2,1,1000]])

K_b = BL4.RigidezGlobal(nn,ne,MAT,ESP,XY,IJ)
#Apenas revisar F(nós)
F_b = BL4.ForcapGlobal(nn,ESP,XY,IJ,nc,P)

K_mod,F_mod = BL4.AplicaCCH(nn,na,AP,K_b,F_b)

U = np.linalg.solve(K_mod,F_mod)
print(U)
Sigma = BL4.CalculaTensaoMalha(nn,ne,MAT,ESP,XY,IJ,U)
