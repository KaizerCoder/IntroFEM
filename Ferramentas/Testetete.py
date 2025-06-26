import numpy as np 

def MontaCEPT(E,v):
    c = E/(1-v**2)
    
    Cv = np.array([[c,c*v,0],[c*v,c,0,],[0,0,c*(1-v)/2]])
    
    return Cv


E = 1e9
v = 0.3
# Aloca Matriz
Ke = np.zeros((8,8))
#Matriz Constitutiva do elemento

Cv = MontaCEPT(E,v)
te = 0.01

#Define os pontos de Gauss

pg = (1/(3**(0.5)))*np.array([[-1,1,1,-1],[-1,-1,1,1]])

ra = (1/(3**(0.5)))*np.array([-1,1,1,-1])
sa = (1/(3**(0.5)))*np.array([-1,-1,1,1])

J = np.array([[0.25,0],[0,0.05]])

#Laço pelos ponos de Gauss
for i in range(4):
    
    r = ra[i]
    s = sa[i]
    
    #Calcula determinante
    dJ = np.linalg.det(J)
        
    #Monta B
    B = np.array([
        [s-1,0,1-s,0,s+1,0,-s-1,0],
        [0,5*r-5,0,-5*r-5,0,5*r+5, 0,5-5*r],
        [5*r-5,s-1,-5*r-5,1-s,5*r+5,s+1,5-5*r,-s-1]
        ])
    
    print(f"B{i}: {B}")

    #Acumula na matriz de rigidez do elemento

    Ke += ((np.transpose(B))@Cv@B)*te*dJ
    

pass
