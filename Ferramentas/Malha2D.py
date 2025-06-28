import gmsh
import meshio
import numpy as np


def Gerar_barra_tracao(filename="barraTracao.msh",L=1.0,h=0.1,nx=5,ny=5,triangulo=False,show_gui=False):
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
    if not triangulo:
        gmsh.model.geo.mesh.setRecombine(2, surf)

    # Physical groups
    gmsh.model.addPhysicalGroup(2, [surf], tag=1)     # MAT1
    gmsh.model.setPhysicalName(2, 1, "MAT1")
    gmsh.model.addPhysicalGroup(1, [l4], tag=2)       # ENGASTE
    gmsh.model.setPhysicalName(1, 2, "Sx")
    gmsh.model.addPhysicalGroup(1, [l1], tag=3)       # Fixação em y da parte inferior
    gmsh.model.setPhysicalName(1, 3, "Sy")
    gmsh.model.addPhysicalGroup(1, [l2], tag=4)       # CARGA
    gmsh.model.setPhysicalName(1, 4, "CARGA")      

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    if show_gui:
        gmsh.fltk.run()
    gmsh.finalize()   

def Gerar_viga_flexao(filename="barraFlexao.msh",L=1.0,h=0.1,nx=10,ny=1,triangulo=False,show_gui=False):
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
    if not triangulo:
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
    if show_gui:
        gmsh.fltk.run()
    gmsh.finalize()

def Gerar_chapa_tracao_Quarto(filename="chapaTracao.msh",L=1.0,R=0.1,p_fino=0.02,p_grosso=0.05,triangulo=False,show_gui=False):
    gmsh.initialize()
    gmsh.model.add("Chapa2D")
    
    #Pontos
    p1 = gmsh.model.geo.addPoint(0,0,0,p_fino)
    p2 = gmsh.model.geo.addPoint(R,0,0,p_fino)
    p3 = gmsh.model.geo.addPoint(L,0,0,p_grosso)
    p4 = gmsh.model.geo.addPoint(L,L,0,p_grosso)
    p5 = gmsh.model.geo.addPoint(0,L,0,p_grosso)
    p6 = gmsh.model.geo.addPoint(0,R,0,p_fino)
    
    #Linhas
    
    l1 = gmsh.model.geo.addLine(p2,p3)
    l2 = gmsh.model.geo.addLine(p3,p4)
    l3 = gmsh.model.geo.addLine(p4,p5)
    l4 = gmsh.model.geo.addLine(p5,p6)
    
    #Raio
    R1  = gmsh.model.geo.add_circle_arc(p6,p1,p2)
    
    #Cria Superficie
    
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4,R1])
    surf = gmsh.model.geo.addPlaneSurface([loop])
    if not triangulo:
        gmsh.model.geo.mesh.setRecombine(2, surf)

    # Physical groups
    gmsh.model.addPhysicalGroup(2, [surf], tag=1)     # MAT1
    gmsh.model.setPhysicalName(2, 1, "MAT1")
    gmsh.model.addPhysicalGroup(1, [l1], tag=2)       # ENGASTE
    gmsh.model.setPhysicalName(1, 2, "Engaste1")
    gmsh.model.addPhysicalGroup(1, [l4], tag=3)       # Fixação em X da parte inferior
    gmsh.model.setPhysicalName(1, 3, "Engaste2")
    gmsh.model.addPhysicalGroup(1, [l2], tag=4)       # CARGA
    gmsh.model.setPhysicalName(1, 4, "CARGA")      
    
    
    # gmsh.mesh
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    if show_gui:
        gmsh.fltk.run()
    gmsh.finalize()  
    
def Gerar_chapa_tracao_completa(nome_arquivo="malha.msh", L=1.0, D=0.4, h_fino=0.02, h_grosso=0.05, triangulo=False, show_gui=False):
    gmsh.initialize()
    gmsh.model.add("chapa_completa")

    R = D / 2
    xc, yc = L / 2, L / 2

    # Furo
    p1 = gmsh.model.occ.addPoint(xc + R, yc, 0, h_fino)
    p2 = gmsh.model.occ.addPoint(xc, yc + R, 0, h_fino)
    p3 = gmsh.model.occ.addPoint(xc - R, yc, 0, h_fino)
    p4 = gmsh.model.occ.addPoint(xc, yc - R, 0, h_fino)
    pc = gmsh.model.occ.addPoint(xc, yc, 0, h_fino)

    c1 = gmsh.model.occ.addCircleArc(p1, pc, p2)
    c2 = gmsh.model.occ.addCircleArc(p2, pc, p3)
    c3 = gmsh.model.occ.addCircleArc(p3, pc, p4)
    c4 = gmsh.model.occ.addCircleArc(p4, pc, p1)

    loop_furo = gmsh.model.occ.addCurveLoop([c1, c2, c3, c4])

    # Chapa
    p5 = gmsh.model.occ.addPoint(0, 0, 0, h_grosso)
    p6 = gmsh.model.occ.addPoint(L, 0, 0, h_grosso)
    p7 = gmsh.model.occ.addPoint(L, L, 0, h_grosso)
    p8 = gmsh.model.occ.addPoint(0, L, 0, h_grosso)

    l1 = gmsh.model.occ.addLine(p5, p6)
    l2 = gmsh.model.occ.addLine(p6, p7)
    l3 = gmsh.model.occ.addLine(p7, p8)
    l4 = gmsh.model.occ.addLine(p8, p5)

    loop_ext = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.occ.addPlaneSurface([loop_ext, loop_furo])

    gmsh.model.occ.synchronize()
        
    gmsh.model.addPhysicalGroup(2, [surf], tag=1)     # MAT1
    gmsh.model.setPhysicalName(2, 1, "MAT1")
    gmsh.model.addPhysicalGroup(1, [l2], tag=2)       # CARGA
    gmsh.model.setPhysicalName(1, 2, "CARGA1")     
    gmsh.model.addPhysicalGroup(1, [l4], tag=3)       # CARGA
    gmsh.model.setPhysicalName(1, 3, "CARGA2")     
    

    # Malha
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", h_fino)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h_grosso)

    if not triangulo:
        gmsh.model.mesh.setRecombine(2, surf)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.mesh.generate(2)
    gmsh.write(nome_arquivo)

    if show_gui:
        gmsh.fltk.run()
    gmsh.finalize()

def AbreVisualizacaoGmsh(filename):
    gmsh.initialize()
    gmsh.open(filename)
    gmsh.fltk.run()
    gmsh.finalize()
    
def Exporta_para_Gmsh(filename,IJ,XY, U=None, TVM=None,detJ=None):
    n_nodes = XY.shape[0]
    n_elements = IJ.shape[0]
    n_nos_por_elem = IJ.shape[1]
    
    if n_nos_por_elem == 3:
        gmsh_elem_type = 2
    elif n_nos_por_elem == 4:
        gmsh_elem_type = 3
    
    node_ids = np.arange(1, n_nodes + 1)
    element_ids = np.arange(1, n_elements + 1)

    with open(filename, "w") as f:
        # Header
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        # Nós
        f.write(f"$Nodes\n{n_nodes}\n")
        for i, (x, y) in enumerate(XY, start=1):
            f.write(f"{i} {x} {y} 0.0\n")
        f.write("$EndNodes\n")

        # Elementos
        f.write(f"$Elements\n{n_elements}\n")
        for eid, conn in zip(element_ids, IJ):
            conn_str = " ".join(str(nid) for nid in conn)
            f.write(f"{eid} {gmsh_elem_type} 2 0 1 {conn_str}\n")
        f.write("$EndElements\n")

        # Deslocamentos 
        if U is not None:
            ux = U[::2]
            uy = U[1::2]
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
        if TVM is not None:    
            f.write("$ElementData\n1\n\"Tensão de von Mises\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_elements))
            for eid, vm in zip(element_ids, TVM):
                f.write(f"{eid} {vm}\n")
            f.write("$EndElementData\n")
        
        # Determinante do Jacobiano
        if detJ is not None:          
            f.write("$ElementData\n1\n\"Determinante do Jacobiano\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_elements))
            for eid, vm in zip(element_ids, detJ):
                f.write(f"{eid} {vm}\n")
            f.write("$EndElementData\n")

def Processar_malha2D(arquivo_msh, materiais, apoios, cargas, espessura):
    mesh = meshio.read(arquivo_msh)

    # 1. Coordenadas dos nós
    XY = mesh.points[:, :2]
    nn = XY.shape[0]

    # 2. Conectividade dos elementos
    if "quad" in mesh.cells_dict:
        elements = mesh.cells_dict["quad"]
        element_type = "quad"
        IJ = elements + 1  # base 1
    elif "triangle" in mesh.cells_dict:
        elements = mesh.cells_dict["triangle"]
        element_type = "triangle"
        IJ = elements + 1  # base 1
    else:
        raise ValueError("A malha deve conter elementos quadriláteros ou triangulares.")

    ne = IJ.shape[0]

    # 3. Materiais
    physical_tags = mesh.cell_data_dict["gmsh:physical"][element_type]
    MAT = []
    ESP = []
    tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}
    for tag in physical_tags:
        nome = tag_to_name.get(tag)
        if nome in materiais:
            E, nu = materiais[nome]
            MAT.append([E, nu])
            ESP.append([espessura])
        else:
            raise ValueError(f"Tag de material '{nome}' não encontrado no dicionário de materiais.")
    MAT = np.array(MAT)
    ESP = np.array(ESP)

    # 4. Apoios
    AP = []
    for name, tipo in apoios.items():
        for entity in mesh.field_data:
            if entity == name:
                tag, dim = mesh.field_data[entity]
                if dim == 1:
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
                    aresta = mesh.cells_dict["line"][i]
                    for e, conn in enumerate(IJ - 1):
                        conn_set = set(conn)
                        if set(aresta).issubset(conn_set):
                            face = None
                            if element_type == "quad":
                                faces = [(0,1),(1,2),(2,3),(3,0)]
                            else:  # triangle
                                faces = [(0,1),(1,2),(2,0)]
                            for f, par in enumerate(faces):
                                if set([conn[par[0]], conn[par[1]]]) == set(aresta):
                                    face = f + 1
                                    break
                            P.append([e + 1, face, direcao, valor])
        elif dim == 0:
            for i, point in enumerate(mesh.cells_dict["vertex"]):
                if mesh.cell_data_dict["gmsh:physical"]["vertex"][i] == tag:
                    if direcao == 1:
                        P.append([point[0] + 1, 1, 0, valor])
                    if direcao == 2:
                        P.append([point[0] + 1, 2, 0, valor])
    P = np.array(P)
    nc = len(P)

    return nn, XY, ne, IJ, MAT, ESP, na, AP, nc, P
