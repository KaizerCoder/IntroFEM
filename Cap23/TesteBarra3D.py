import numpy as np
import Elementos3D as H3D
import sympy as sp
import meshio
import Malha3D


#Entrada de Dados 

# #Gera Malha(básica)
arquivo = 'Barra3D.msh'

Materiais = {"MAT1":[1e9,0.3]}
apoios = {'ENGASTE':'xyz','ENGASTE2':"yz"}
cargas = [("CARGA",1,1000)]


Malha3D.Gerar_barra_tracao3D(arquivo,L=1,h=0.1,t=0.1,nx=10,ny=10,nz=10,show_gui=False)
nn, XYZ, ne, IJ, MAT, na, AP, nc, P = Malha3D.Processar_malha3D(arquivo,Materiais,apoios,cargas)

#Monta Rigidez Global
K_b = H3D.RigidezGlobal(nn,ne,MAT,XYZ,IJ)
#Monta Força Global
F_b = H3D.ForcapGlobal(nn,XYZ,IJ,nc,P)
#Aplica as CChs homogêneas
K_mod,F_mod = H3D.AplicaCCH(nn,na,AP,K_b,F_b)

#Calcula o Deslocamento
U = np.linalg.solve(K_mod,F_mod)

#Calcula a Tensão ao longo dos elementos
Sigma = H3D.CalculaTensaoMalha(nn,ne,MAT,XYZ,IJ,U)

# #Calcula Von Mises
TVM = H3D.TensaoVonMises3D(ne,Sigma)

Malha3D.Exporta_para_Gmsh(arquivo,IJ, XYZ, U,TVM)
Malha3D.AbreVisualizacaoGmsh(arquivo)
# # Malha2D.Exporta_e_abre_gmsh(arquivo,IJ, XY, U, TVM,detJ)
