import numpy as np
import Elementos2D as BL4
import sympy as sp
import meshio
import Malha2D


#Entrada de Dados 

# #Gera Malha(básica)
arquivo = 'ChapaCompleta.msh'

Materiais = {"MAT1":[1e9,0.3]}
apoios = {}
cargas = [("CARGA1",1,1000),("CARGA2",1,1000)]


Malha2D.Gerar_chapa_tracao_completa(nome_arquivo=arquivo, L=1.0, D=0.1, h_fino=0.01, h_grosso=0.05)
nn, XY, ne, IJ, MAT, ESP, na, AP, nc, P = Malha2D.Processar_malha2D(arquivo,Materiais,apoios,cargas,espessura=0.01)

# #Entrada de Dados 

# arquivo = 'Triangulo.msh'
# #Número de nós
# nn = 6
# XY = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.0],[0.0,0.1],[0.5,0.1],[1.0,0.1]])
# ne = 4
# IJ = np.array([[1,2,5],[2,3,6],[1,5,4],[2,6,5]])

# MAT = np.array([[1e9,0.3],[1e9,0.3],[1e9,0.3],[1e9,0.3]])

# ESP = np.array([[0.01],[0.01],[0.01],[0.01]])

# na = 5
# AP = np.array([[1,1,0],[1,2,0],[2,2,0],[3,2,0],[4,1,0]])

# nc=1
# P = np.array([[2,2,1,1000]])



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

detJ = BL4.CalculaDetJ(ne,XY,IJ)

#Calcula Von Mises
TVM = BL4.TensaoVonMises2D(ne,Sigma)

Malha2D.Exporta_para_Gmsh(arquivo,IJ, XY, U, TVM,detJ)
Malha2D.AbreVisualizacaoGmsh(arquivo)
# Malha2D.Exporta_e_abre_gmsh(arquivo,IJ, XY, U, TVM,detJ)
