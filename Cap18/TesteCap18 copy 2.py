import numpy as np
import Elementos2D as BL4
import sympy as sp
import meshio
import Malha2D


#Entrada de Dados 

#Gera Malha(básica)




#Entrada de Dados 

Materiais = {"MAT1":[10e3,0]}
apoios = {"ENGASTE":"engaste"}
cargas = [("CARGA",2,1)]

n = 2
arquivo = "barraFlexaoBolha.msh"
Malha2D.Gerar_viga_flexao(arquivo,L=1,h=0.1,nx=10*n,ny=n)
nn, XY, ne, IJ, MAT, ESP, na, AP, nc, P = Malha2D.Processar_malha2D(arquivo,Materiais,apoios,cargas,espessura=0.1)

#Monta Rigidez Global
K_b = BL4.RigidezGlobal(nn,ne,MAT,ESP,XY,IJ,bolha=True)
#Monta Força Global
F_b = BL4.ForcapGlobal(nn,ESP,XY,IJ,nc,P)
#Aplica as CChs homogêneas
K_mod,F_mod = BL4.AplicaCCH(nn,na,AP,K_b,F_b)

#Calcula o Deslocamento
U = np.linalg.solve(K_mod,F_mod)
print(U[5])
# sp.pprint(sp.latex(sp.Matrix(U)))


#Calcula a Tensão ao longo dos elementos
Sigma = BL4.CalculaTensaoMalha(nn,ne,MAT,ESP,XY,IJ,U,bolha=True)


detJ = BL4.CalculaDetJ(ne,XY,IJ)

#Calcula Von Mises
TVM = BL4.TensaoVonMises2D(ne,Sigma)

Malha2D.export_to_gmsh_post(arquivo, XY, U, TVM,IJ,detJ)