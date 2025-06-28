import numpy as np
import Elementos2D as BL4
import sympy as sp
import meshio
import Malha2D


#Entrada de Dados 

# #Gera Malha(básica)
Materiais = {"MAT1":[1e9,0.3]}
apoios = {"Engaste1":"engaste","Engaste2":"engaste"}
cargas = [("CARGA",1,1000)]
arquivo = "ChapaTracao.msh"
Malha2D.gerar_chapa_tracao("ChapaTracao.msh",L=1.0,h=0.1,R=0.3,p_fino=0.01,p_grosso=0.05)
nn, XY, ne, IJ, MAT, ESP, na, AP, nc, P = Malha2D.processar_malha2D(arquivo,Materiais,apoios,cargas,espessura=0.01)

#Monta Rigidez Global
K_b = BL4.RigidezGlobal(nn,ne,MAT,ESP,XY,IJ)
#Monta Força Global
F_b = BL4.ForcapGlobal(nn,ESP,XY,IJ,nc,P)
#Aplica as CChs homogêneas
K_mod,F_mod = BL4.AplicaCCH(nn,na,AP,K_b,F_b)

#Calcula o Deslocamento
U = np.linalg.solve(K_mod,F_mod)
# print(U)
sp.pprint(sp.latex(sp.Matrix(U)))

UX = U[0::2]
UY = U[1::2]

detJ = BL4.CalculaDetJ(ne,XY,IJ)

#Calcula a Tensão ao longo dos elementos
Sigma = BL4.CalculaTensaoMalha(nn,ne,MAT,ESP,XY,IJ,U)


#Calcula Von Mises
TVM = BL4.TensaoVonMises2D(ne,Sigma)


Malha2D.export_to_gmsh_post(arquivo, XY, U, TVM,IJ,detJ)
