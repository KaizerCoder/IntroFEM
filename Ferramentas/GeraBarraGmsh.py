import gmsh
import meshio
import numpy as np


def gerar_barra_tracao(filename="barraTracao.msh",L=1.0,h=0.1,nx=5,ny=5):
    gmsh.initialize()
    gmsh.model.add("Barra2D")
    
    #Pontos
    p1 = gmsh.model.geo.addPoint(0,0,0)
    p2 = gmsh.model.geo.addPoint(L,0,0)
    p3 = gmsh.model.geo.addPoint(L,h,0)
    p4 = gmsh.model.geo.addPoint(0,h,0)
    
    #Linhas
    
    l1 = gmsh.model.geo.addLine(p1,p2)
    l2 = gmsh.model.geo.addLine(p2,p3)
    l3 = gmsh.model.geo.addLine(p3,p4)
    l4 = gmsh.model.geo.addLine(p4,p1)
    
    #Cria Superficie
    
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])
    
    # Transfinite (malha estruturada)
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, nx + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, nx + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, ny + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, ny + 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(surf)
    gmsh.model.geo.mesh.setRecombine(2, surf)

    # Physical groups
    gmsh.model.addPhysicalGroup(2, [surf], tag=1)     # MAT1
    gmsh.model.setPhysicalName(2, 1, "MAT1")
    gmsh.model.addPhysicalGroup(1, [l4], tag=2)       # ENGASTE
    gmsh.model.setPhysicalName(1, 2, "Sx")
    gmsh.model.addPhysicalGroup(1, [l1], tag=3)       # Fixação em X da parte inferior
    gmsh.model.setPhysicalName(1, 3, "Sy")
    gmsh.model.addPhysicalGroup(1, [l2], tag=4)       # CARGA
    gmsh.model.setPhysicalName(1, 4, "CARGA")      

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    gmsh.finalize()
    

def gerar_barra_api_flexao(filename="barraFlexao.msh",L=1.0,h=0.1,nx=10,ny=1):
    gmsh.initialize()
    gmsh.model.add("Barra2D")
    
    #Pontos
    p1 = gmsh.model.geo.addPoint(0,0,0)
    p2 = gmsh.model.geo.addPoint(L,0,0)
    p3 = gmsh.model.geo.addPoint(L,h,0)
    p4 = gmsh.model.geo.addPoint(0,h,0)
    
    #Linhas
    
    l1 = gmsh.model.geo.addLine(p1,p2)
    l2 = gmsh.model.geo.addLine(p2,p3)
    l3 = gmsh.model.geo.addLine(p3,p4)
    l4 = gmsh.model.geo.addLine(p4,p1)
    
    #Cria Superficie
    
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])
    
    # Transfinite (malha estruturada)
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, nx + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, nx + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, ny + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, ny + 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(surf)
    gmsh.model.geo.mesh.setRecombine(2, surf)

    # Physical groups
    gmsh.model.addPhysicalGroup(2, [surf], tag=1)     # MAT1
    gmsh.model.setPhysicalName(2, 1, "MAT1")
    gmsh.model.addPhysicalGroup(1, [l4], tag=2)       # ENGASTE
    gmsh.model.setPhysicalName(1, 2, "ENGASTE")
    gmsh.model.addPhysicalGroup(0, [p3], tag=3)       # CARGA
    gmsh.model.setPhysicalName(0, 3, "CARGA")      

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    gmsh.finalize()
    
        
def export_to_gmsh_post(filename, XY, U, TVM,IJ):
    n_nodes = XY.shape[0]
    n_elements = IJ.shape[0]
    
    node_ids = np.arange(1, n_nodes + 1)
    element_ids = np.arange(1, n_elements + 1)

    ux = U[::2]
    uy = U[1::2]

    with open(filename, "w") as f:
        # Mesh header
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        # Nodes
        f.write(f"$Nodes\n{n_nodes}\n")
        for i, (x, y) in enumerate(XY, start=1):
            f.write(f"{i} {x} {y} 0.0\n")
        f.write("$EndNodes\n")

        # Elements (quadrilateral 4-node)
        f.write(f"$Elements\n{n_elements}\n")
        for eid, conn in zip(element_ids, IJ):
            conn_str = " ".join(str(nid) for nid in conn)
            f.write(f"{eid} 3 2 0 1 {conn_str}\n")  # type 3 = quad, 2 tags
        f.write("$EndElements\n")

        # Deslocamento Ux
        f.write("$NodeData\n1\n\"Deslocamento Ux\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
        for nid, val in zip(node_ids, ux):
            f.write(f"{int(nid)} {float(val)}\n")
        f.write("$EndNodeData\n")

        # Deslocamento Uy
        f.write("$NodeData\n1\n\"Deslocamento Uy\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
        for nid, val in zip(node_ids, uy):
            f.write(f"{int(nid)} {float(val)}\n")
        f.write("$EndNodeData\n")

        # Tensão de Von Mises
        f.write("$ElementData\n1\n\"Tensão de von Mises\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_elements))
        for eid, vm in zip(element_ids, TVM):
            f.write(f"{eid} {vm}\n")
        f.write("$EndElementData\n")

def processar_malha2D(arquivo_msh, materiais, apoios, cargas,espessura):
    mesh = meshio.read(arquivo_msh)

    # 1. Coordenadas dos nós
    XY = mesh.points[:, :2]
    nn = XY.shape[0]

    # 2. Conectividade dos elementos quadriláteros
    elements = mesh.cells_dict.get("quad", None)
    if elements is None:
        raise ValueError("A malha deve conter elementos quadriláteros.")
    IJ = elements + 1  # base 1
    ne = IJ.shape[0]

    # 3. Materiais - Associa cada elemento ao seu material com base no physical tag
    physical_tags = mesh.cell_data_dict["gmsh:physical"]["quad"]
    MAT = []
    ESP = []
    tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}  # tag -> nome
    for tag in physical_tags:
        nome = tag_to_name.get(tag)
        if nome in materiais:
            E, nu = materiais[nome]
            MAT.append([E, nu])
            ESP.append([espessura])  # Aqui você pode colocar espessura variável por material se desejar
        else:
            raise ValueError(f"Tag de material '{nome}' não encontrado no dicionário de materiais.")
    MAT = np.array(MAT)
    ESP = np.array(ESP)

    # 4. Apoios: associar nodes com restrições
    AP = []
    for name, tipo in apoios.items():
        for entity in mesh.field_data:
            if entity == name:
                tag, dim = mesh.field_data[entity]
                if dim == 1:
                    lines = []
                    lines = mesh.cells_dict.get("line", [])
                    phys_tags = mesh.cell_data_dict["gmsh:physical"]["line"]
                    for i, t in enumerate(phys_tags):
                        if t == tag:
                            for node in lines[i]:
                                if tipo in ("engaste", "x"):
                                    AP.append([node + 1, 1, 0])
                                if tipo in ("engaste", "y"):
                                    AP.append([node + 1, 2, 0])
                elif dim == 0:
                    # ponto individual
                    for i, point in enumerate(mesh.cells_dict["vertex"]):
                        if mesh.cell_data_dict["gmsh:physical"]["vertex"][i] == tag:
                            if tipo in ("x", "engaste"):
                                AP.append([point[0] + 1, 1, 0])
                            if tipo in ("y", "engaste"):
                                AP.append([point[0] + 1, 2, 0])
    AP = np.array(AP)
    na = len(AP)

    # 5. Cargas
    P = []
    for nome, direcao, valor in cargas:
        tag, dim = mesh.field_data[nome]
        if dim == 1:
            for i, t in enumerate(mesh.cell_data_dict["gmsh:physical"]["line"]):
                if t == tag:
                    # associar ao elemento por proximidade dos nós da aresta
                    aresta = mesh.cells_dict["line"][i]
                    for e, conn in enumerate(IJ - 1):
                        if set(aresta).issubset(set(conn)):
                            # identificar qual face é (1 a 4)
                            face = None
                            for f, par in enumerate([(0,1),(1,2),(2,3),(3,0)]):
                                if set([conn[par[0]], conn[par[1]]]) == set(aresta):
                                    face = f + 1
                                    break
                            P.append([e + 1, face, direcao, valor])
        elif dim == 0:
            # ponto individual
            for i, point in enumerate(mesh.cells_dict["vertex"]):
                if mesh.cell_data_dict["gmsh:physical"]["vertex"][i] == tag:
                    if direcao==1:
                        P.append([point[0] + 1,1,0, valor])
                    if direcao==2:
                        P.append([point[0] + 1,2,0, valor])
    P = np.array(P)
    nc = len(P)

    return nn, XY, ne, IJ, MAT, ESP, na, AP, nc, P