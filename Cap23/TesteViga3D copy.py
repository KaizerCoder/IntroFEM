import numpy as np
import Elementos3D as H3D
import sympy as sp
import meshio
import Malha3D


#Entrada de Dados 

# #Gera Malha(básica)
arquivo = 'Viga3D.msh'

Materiais = {"MAT1":[10e3,0]}
apoios = {'ENGASTE':'xyz'}
cargas = [("CARGA",2,100),("CARGA",3,100)]

n = 1
Malha3D.Gerar_Viga_Flexao3D_Corrigido(arquivo,L=1,h=0.1,t=0.1,nx=10*n,ny=n,nz=n,show_gui=False)
nn, XYZ, ne, IJ, MAT, na, AP, nc, P = Malha3D.Processar_malha3D(arquivo,Materiais,apoios,cargas)

#Monta Rigidez Global
K_b = H3D.RigidezGlobal(nn,ne,MAT,XYZ,IJ,bolha=True)
#Monta Força Global
F_b = H3D.ForcapGlobal(nn,XYZ,IJ,nc,P)
#Aplica as CChs homogêneas
K_mod,F_mod = H3D.AplicaCCH(nn,na,AP,K_b,F_b)

#Calcula o Deslocamento
U = np.linalg.solve(K_mod,F_mod)

#Calcula a Tensão ao longo dos elementos
Sigma = H3D.CalculaTensaoMalha(nn,ne,MAT,XYZ,IJ,U,bolha=True)

# #Calcula Von Mises
TVM = H3D.TensaoVonMises3D(ne,Sigma)

Malha3D.Exporta_para_Gmsh_GEO(arquivo,IJ, XYZ, U,TVM,None,F_b,P)
Malha3D.AbreVisualizacaoGmsh(arquivo)
