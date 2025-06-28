import numpy as np
import Elementos2D as BL4
import sympy as sp
import meshio
import Malha2D


#Entrada de Dados 

#Gera Malha(básica)
Materiais = {"MAT1":[1e9,0.3]}
apoios = {"Sx":"x","Sy":"y"}
cargas = [("CARGA",1,1000)]

arquivo = "barraTracao3.msh"

Malha2D.Gerar_barra_tracao(arquivo,L=1,h=0.1,nx=2,ny=2)
nn, XY, ne, IJ, MAT, ESP, na, AP, nc, P = Malha2D.Processar_malha2D(arquivo,Materiais,apoios,cargas,espessura=0.1)

#Monta Rigidez Global
K_b = BL4.RigidezGlobal(nn,ne,MAT,ESP,XY,IJ)
#Monta Força Global
F_b = BL4.ForcapGlobal(nn,ESP,XY,IJ,nc,P)
#Aplica as CChs homogêneas
K_mod,F_mod = BL4.AplicaCCH(nn,na,AP,K_b,F_b)

#Calcula o Deslocamento
U = np.linalg.solve(K_mod,F_mod)
print(U)
sp.pprint(sp.latex(sp.Matrix(U)))

UX = U[0::2]
UY = U[1::2]


#Calcula a Tensão ao longo dos elementos
Sigma = BL4.CalculaTensaoMalha(nn,ne,MAT,ESP,XY,IJ,U)

#Calcula Von Mises
TVM = BL4.TensaoVonMises2D(ne,Sigma)


Malha2D.Exporta_para_Gmsh(arquivo,IJ, XY, U, TVM)

Malha2D.AbreVisualizacaoGmsh(arquivo)