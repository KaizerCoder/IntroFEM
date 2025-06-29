import gmsh
import meshio
import numpy as np


# --- Obter centroide de cada face manualmente
def get_face_centroid(face_tag):
    types, _, node_tags = gmsh.model.mesh.getElements(2, face_tag)
    node_ids = node_tags[0]
    all_node_ids, coords, _ = gmsh.model.mesh.getNodes()
    coord_dict = {nid: coords[i*3:i*3+3] for i, nid in enumerate(all_node_ids)}
    coords_face = [coord_dict[nid] for nid in node_ids if nid in coord_dict]
    cx = sum(pt[0] for pt in coords_face) / len(coords_face)
    cy = sum(pt[1] for pt in coords_face) / len(coords_face)
    cz = sum(pt[2] for pt in coords_face) / len(coords_face)
    return (cx, cy, cz)

def Gerar_barra_tracao3D(filename="barraTracao3D.msh",L=1.0,h=0.1,t=0.1,nx=5,ny=5,nz=5,show_gui=False):
    gmsh.initialize()
    gmsh.model.add("Barra3D")
    
    #Pontos do plano XY
    p1 = gmsh.model.geo.addPoint(0,0,0)
    p2 = gmsh.model.geo.addPoint(L,0,0)
    p3 = gmsh.model.geo.addPoint(L,h,0)
    p4 = gmsh.model.geo.addPoint(0,h,0)
    
    #Linhas do plano XY
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
    
    # --- Extrusão em Z (espessura)
    extrude_entities = gmsh.model.geo.extrude([(2, surf)], 0, 0, t, numElements=[nz], recombine=True)

    # --- Obter volume (primeiro elemento 3D retornado)
    volume = [e[1] for e in extrude_entities if e[0] == 3][0]
    gmsh.model.addPhysicalGroup(3, [volume], tag=1)
    gmsh.model.setPhysicalName(3, 1, "MAT1")

    # --- Sincronizar antes de trabalhar com malha
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    # --- Identificar faces laterais (dimension 2)
    face_tags = [e[1] for e in extrude_entities if e[0] == 2]

    # --- Criar Physical Groups nas faces
    for face in face_tags:
        cx, cy, cz = get_face_centroid(face)
        if abs(cx) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [face], tag=2)
            gmsh.model.setPhysicalName(2, 2, "ENGASTE")
            break
    

    for face in face_tags:
        cx, cy, cz = get_face_centroid(face)
        if abs(cx - L) < 1e-6:
            gmsh.model.addPhysicalGroup(2, [face], tag=3)
            gmsh.model.setPhysicalName(2, 3, "CARGA")
            break
    
    
    for face in face_tags:
        cx, cy, cz = get_face_centroid(face)
        if abs(cx-L/2) < 1e-6 and (cy)<1e-6 and (cz-t/2)<1e6-6:
            gmsh.model.addPhysicalGroup(2, [face], tag=4)
            gmsh.model.setPhysicalName(2, 4, "ENGASTE2")
            break    
    
    gmsh.model.geo.synchronize()

    gmsh.write(filename)
    if show_gui:
        gmsh.fltk.run()
    gmsh.finalize()   

# def Gerar_Viga_Flexao3D(filename="barraTracao3D.msh",L=1.0,h=0.1,t=0.1,nx=5,ny=5,nz=5,show_gui=False):
#     gmsh.initialize()
#     gmsh.model.add("Viga3D")
    
#     #Pontos do plano XY
#     p1 = gmsh.model.geo.addPoint(0,0,0)
#     p2 = gmsh.model.geo.addPoint(L,0,0)
#     p3 = gmsh.model.geo.addPoint(L,h,0)
#     p4 = gmsh.model.geo.addPoint(0,h,0)
    
#     #Linhas do plano XY
#     l1 = gmsh.model.geo.addLine(p1,p2)
#     l2 = gmsh.model.geo.addLine(p2,p3)
#     l3 = gmsh.model.geo.addLine(p3,p4)
#     l4 = gmsh.model.geo.addLine(p4,p1)
    
#     #Cria Superficie
    
#     loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
#     surf = gmsh.model.geo.addPlaneSurface([loop])
    
#     # Transfinite (malha estruturada)
#     gmsh.model.geo.mesh.setTransfiniteCurve(l1, nx + 1)
#     gmsh.model.geo.mesh.setTransfiniteCurve(l3, nx + 1)
#     gmsh.model.geo.mesh.setTransfiniteCurve(l2, ny + 1)
#     gmsh.model.geo.mesh.setTransfiniteCurve(l4, ny + 1)
#     gmsh.model.geo.mesh.setTransfiniteSurface(surf)
#     gmsh.model.geo.mesh.setRecombine(2, surf)
    
#     # --- Extrusão em Z (espessura)
#     extrude_entities = gmsh.model.geo.extrude([(2, surf)], 0, 0, t, numElements=[nz], recombine=True)

#     # --- Obter volume (primeiro elemento 3D retornado)
#     volume = [e[1] for e in extrude_entities if e[0] == 3][0]
#     gmsh.model.addPhysicalGroup(3, [volume], tag=1)
#     gmsh.model.setPhysicalName(3, 1, "MAT1")

#     # --- Sincronizar antes de trabalhar com malha
#     gmsh.model.geo.synchronize()
#     gmsh.model.mesh.generate(3)

#     # --- Identificar faces laterais (dimension 2)
#     face_tags = [e[1] for e in extrude_entities if e[0] == 2]

#     # --- Criar Physical Groups nas faces
#     for face in face_tags:
#         cx, cy, cz = get_face_centroid(face)
#         if abs(cx) < 1e-6:
#             gmsh.model.addPhysicalGroup(2, [face], tag=2)
#             gmsh.model.setPhysicalName(2, 2, "ENGASTE")
#             break
#         elif abs(cx - L) < 1e-6:
#             carga_face = face
    

#     for face in face_tags:
#         cx, cy, cz = get_face_centroid(face)
#         if abs(cx - L) < 1e-6:
#             gmsh.model.addPhysicalGroup(2, [face], tag=3)
#             gmsh.model.setPhysicalName(2, 3, "CARGA")
#             break
    
    
#     for face in face_tags:
#         cx, cy, cz = get_face_centroid(face)
#         if abs(cx-L/2) < 1e-6 and (cy)<1e-6 and (cz-t/2)<1e6-6:
#             gmsh.model.addPhysicalGroup(2, [face], tag=4)
#             gmsh.model.setPhysicalName(2, 4, "ENGASTE2")
#             break    
    
#     gmsh.model.geo.synchronize()

#     gmsh.write(filename)
#     if show_gui:
#         gmsh.fltk.run()
#     gmsh.finalize()   


# def Gerar_Viga_Flexao3D2(filename="vigaFlexao3D.msh",L=1.0,h=0.1,t=0.1,nx=5,ny=5,nz=5,show_gui=False):
#     # Criar retângulo base no plano XY (z=0)
#     p1 = gmsh.model.geo.addPoint(0, 0, 0)
#     p2 = gmsh.model.geo.addPoint(L, 0, 0)
#     p3 = gmsh.model.geo.addPoint(L, h, 0)
#     p4 = gmsh.model.geo.addPoint(0, h, 0)

#     l1 = gmsh.model.geo.addLine(p1, p2)
#     l2 = gmsh.model.geo.addLine(p2, p3)
#     l3 = gmsh.model.geo.addLine(p3, p4)
#     l4 = gmsh.model.geo.addLine(p4, p1)
    
#     cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
#     surf = gmsh.model.geo.addPlaneSurface([cl])
    
#     # X-direção (p1–p2 e p4–p3)
#     gmsh.model.geo.mesh.setTransfiniteCurve(l1, nx + 1)  # p1 → p2
#     gmsh.model.geo.mesh.setTransfiniteCurve(l3, nx + 1)  # p3 → p4

#     # Y-direção (p2–p3 e p1–p4)
#     gmsh.model.geo.mesh.setTransfiniteCurve(l2, ny + 1)  # p2 → p3
#     gmsh.model.geo.mesh.setTransfiniteCurve(l4, ny + 1)  # p4 → p1

#     # Superfície 2D estruturada
#     gmsh.model.geo.mesh.setTransfiniteSurface(surf)
#     gmsh.model.geo.mesh.setRecombine(2, surf)  # para hexaedros depois da extrusão
    
    
#     gmsh.model.geo.synchronize()

#     # Extrude para cima (z positivo)
#     entities_up = gmsh.model.geo.extrude([(2, surf)], 0, 0, t/2, numElements=[nz], recombine=True)
#     vol_up = [e[1] for e in entities_up if e[0] == 3][0]

#     # Extrude para baixo (z negativo)
#     entities_down = gmsh.model.geo.extrude([(2, surf)], 0, 0, -t/2, numElements=[nz], recombine=True)
#     vol_down = [e[1] for e in entities_down if e[0] == 3][0]
    
#     gmsh.model.geo.synchronize()

#     # União booleana dos volumes up e down
#     result = gmsh.model.occ.fuse([(3, vol_up)], [(3, vol_down)])
#     gmsh.model.geo.synchronize()

#     vol = result[0][1]

#     # Criar Physical Group do volume - material
#     gmsh.model.addPhysicalGroup(3, [vol], tag=1)
#     gmsh.model.setPhysicalName(3, 1, "MAT1")

#     # --- Identificar faces laterais (dimension 2)
#     face_tags = [e[1] for e in vol if e[0] == 2]

#     # --- Criar Physical Groups nas faces
#     for face in face_tags:
#         cx, cy, cz = get_face_centroid(face)
#         if abs(cx) < 1e-6:
#             gmsh.model.addPhysicalGroup(2, [face], tag=2)
#             gmsh.model.setPhysicalName(2, 2, "ENGASTE")
#             break

#     gmsh.model.addPhysicalGroup(0, [p3], tag=3)
#     gmsh.model.setPhysicalName(0, 3, "CARGA")
    
    
#     gmsh.write(filename)
#     gmsh.finalize()


def Gerar_Viga_Flexao3D_Corrigido(
    filename="viga_3d_corrigida.msh", 
    L=1.0, 
    h=0.1, 
    t=0.1, 
    nx=20, 
    ny=4, 
    nz=4, 
    show_gui=False
):
    """
    Gera uma malha 3D de uma viga engastada em uma extremidade com carga pontual na outra.
    Versão corrigida e simplificada.

    Args:
        filename (str): Nome do arquivo .msh de saída.
        L (float): Comprimento da viga (direção x).
        h (float): Altura da viga (direção y).
        t (float): Espessura da viga (direção z).
        nx (int): Número de elementos ao longo do comprimento.
        ny (int): Número de elementos ao longo da altura.
        nz (int): Número de elementos ao longo da espessura.
        show_gui (bool): Se True, exibe a GUI do Gmsh ao final.
    """
    gmsh.initialize()
    gmsh.model.add("viga_3d")

    # --- 1. Definição da Geometria (Base 2D) ---
    # Para centralizar a viga no eixo z=0, criamos a base no plano z = -t/2.
    z_base = -t / 2
    p1 = gmsh.model.geo.addPoint(0, 0, z_base)
    p2 = gmsh.model.geo.addPoint(L, 0, z_base)
    p3 = gmsh.model.geo.addPoint(L, h, z_base)
    p4 = gmsh.model.geo.addPoint(0, h, z_base)

    # Linhas da face base
    l_inf = gmsh.model.geo.addLine(p1, p2) # Linha inferior
    l_dir = gmsh.model.geo.addLine(p2, p3) # Linha direita (x=L)
    l_sup = gmsh.model.geo.addLine(p3, p4) # Linha superior
    l_esq = gmsh.model.geo.addLine(p4, p1) # Linha esquerda (x=0, onde será o engaste)

    # Loop e Superfície
    cl = gmsh.model.geo.addCurveLoop([l_inf, l_dir, l_sup, l_esq])
    superficie_base = gmsh.model.geo.addPlaneSurface([cl])

    # --- 2. Malha Estruturada (Transfinite) ---
    # Definimos o número de nós nas arestas da superfície base.
    # Isso controla a densidade da malha em cada direção.
    gmsh.model.geo.mesh.setTransfiniteCurve(l_inf, nx + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_sup, nx + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_dir, ny + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_esq, ny + 1)

    # Marcamos a superfície como transfinita e pedimos para recombinar os
    # triângulos em quadriláteros, que se tornarão hexaedros após a extrusão.
    gmsh.model.geo.mesh.setTransfiniteSurface(superficie_base)
    gmsh.model.geo.mesh.setRecombine(2, superficie_base)
    
    # Sincronizamos o modelo geométrico antes da extrusão
    gmsh.model.geo.synchronize()

    # --- 3. Extrusão para criar o Volume 3D ---
    # Realizamos UMA ÚNICA extrusão da superfície base ao longo da espessura 't'.
    # A função extrude retorna uma lista com as novas entidades criadas:
    # O primeiro item é o volume, o segundo é a face do topo, e os seguintes
    # são as faces laterais, na mesma ordem das linhas do loop original.
    entidades_extrudadas = gmsh.model.geo.extrude(
        [(2, superficie_base)], 0, 0, t, numElements=[nz], recombine=True
    )
    
    # Sincronizamos novamente após a operação principal
    gmsh.model.geo.synchronize()

    # --- 4. Definição dos Grupos Físicos (Physical Groups) ---
    # Grupos Físicos são essenciais para aplicar materiais e condições de contorno
    # em softwares de Elementos Finitos (FEA).

    # Grupo Físico para o Volume (para atribuição do material)
    tag_volume = [e[1] for e in entidades_extrudadas if e[0] == 3][0]
    gmsh.model.addPhysicalGroup(3, [tag_volume], tag=1)
    gmsh.model.setPhysicalName(3, 1, "MAT1")

    # Grupo Físico para a Face de Engaste (em x=0)
    # A face do engaste é criada pela extrusão da linha l_esq.
    # Pela ordem de retorno de 'extrude', sabemos qual é a tag dessa face.
    # [0]=Volume, [1]=Face Topo, [2]=Face de l_inf, [3]=Face de l_dir, [4]=Face de l_sup, [5]=Face de l_esq
    tag_face_engaste = entidades_extrudadas[5][1]
    gmsh.model.addPhysicalGroup(2, [tag_face_engaste], tag=2)
    gmsh.model.setPhysicalName(2, 2, "ENGASTE")

    # Grupo Físico para o Ponto de Carga
    # Para aplicar uma carga pontual em (L, h, 0), que fica no meio da aresta superior,
    # criamos um ponto explicitamente e o forçamos a fazer parte da malha do volume.
    ponto_carga = gmsh.model.geo.addPoint(L, h, 0)
    gmsh.model.geo.synchronize() # Sincroniza para o novo ponto ser reconhecido
    # "Embed" força o gerador de malha a criar um nó na localização deste ponto.
    # gmsh.model.mesh.embed(0, [ponto_carga], 3, 1)
    
    tag_face_carga = entidades_extrudadas[3][1]
    gmsh.model.addPhysicalGroup(2, [tag_face_carga], tag=3)
    gmsh.model.setPhysicalName(2, 3, "CARGA")
    
    # for face in entidades_extrudadas:
    #     cx, cy, cz = get_face_centroid(face[1])
    #     if abs(cx - L) < 1e-6:
    #         gmsh.model.addPhysicalGroup(2, [face[1]], tag=3)
    #         gmsh.model.setPhysicalName(2, 3, "CARGA")
    #         break
    
    
    # --- 5. Geração da Malha e Salvamento ---
    gmsh.model.mesh.generate(3) # Gera a malha 3D
    gmsh.write(filename) # Salva o arquivo .msh

    # Abre a interface gráfica do Gmsh, se solicitado
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
    
def Exporta_para_Gmsh(filename,IJ,XY, U=None, TVM=None,detJ=None,F=None,P = None,AP=None):
    n_nodes = XY.shape[0]
    n_elements = IJ.shape[0]
    n_nos_por_elem = IJ.shape[1]
    
    if n_nos_por_elem == 3:
        gmsh_elem_type = 2
    elif n_nos_por_elem == 4:
        gmsh_elem_type = 3
    elif n_nos_por_elem == 8:
        gmsh_elem_type = 5
        
    node_ids = np.arange(1, n_nodes + 1)
    element_ids = np.arange(1, n_elements + 1)

    with open(filename, "w") as f:
        # Header
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        # Nós
        f.write(f"$Nodes\n{n_nodes}\n")
        for i, (x, y, z) in enumerate(XY, start=1):
            f.write(f"{i} {x} {y} {z}\n")
        f.write("$EndNodes\n")

        # Elementos
        f.write(f"$Elements\n{n_elements}\n")
        for eid, conn in zip(element_ids, IJ):
            conn_str = " ".join(str(nid) for nid in conn)
            f.write(f"{eid} {gmsh_elem_type} 2 0 1 {conn_str}\n")
        f.write("$EndElements\n")

        # Deslocamentos
        if U is not None:
            ux = U[0::3]
            uy = U[1::3]
            uz = U[2::3]

            f.write("$NodeData\n1\n\"Deslocamento Ux\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
            for nid, val in zip(node_ids, ux):
                f.write(f"{int(nid)} {float(val)}\n")
            f.write("$EndNodeData\n")

            f.write("$NodeData\n1\n\"Deslocamento Uy\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
            for nid, val in zip(node_ids, uy):
                f.write(f"{int(nid)} {float(val)}\n")
            f.write("$EndNodeData\n")

            f.write("$NodeData\n1\n\"Deslocamento Uz\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
            for nid, val in zip(node_ids, uz):
                f.write(f"{int(nid)} {float(val)}\n")
            f.write("$EndNodeData\n")

         # Vetor de Forças Nodais
        if F is not None:
            fx = F[0::3]*100
            fy = F[1::3]*100
            fz = F[2::3]*100

            f.write("$NodeData\n1\n\"Vetores de Força\"\n1\n0.0\n3\n0\n3\n{}\n".format(n_nodes))
            for nid, x, y, z in zip(node_ids, fx, fy, fz):
                f.write(f"{nid} {x} {y} {z}\n")
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
            for eid, val in zip(element_ids, detJ):
                f.write(f"{eid} {val}\n")
            f.write("$EndElementData\n")
            
        if P is not None and IJ.shape[1] == 8:
            force_field = np.zeros((XY.shape[0], 3))  # (n_nodes, [fx, fy, fz])

            face_map = [
                (0, 1, 2, 3),  # Face 1: bottom
                (4, 5, 6, 7),  # Face 2: top
                (0, 1, 5, 4),  # Face 3: front
                (1, 2, 6, 5),  # Face 4: right
                (2, 3, 7, 6),  # Face 5: back
                (3, 0, 4, 7)   # Face 6: left
            ]

            for eid, face_id, direcao, valor in P:
                conn = IJ[int(eid) - 1] - 1  # base 0
                face_nodes = [conn[i] for i in face_map[int(face_id) - 1]]
                for nid in face_nodes:
                    force_field[nid, direcao - 1] += (valor * 1) / len(face_nodes)

            f.write("$NodeData\n1\n\"Carga Aplicada\"\n1\n0.0\n3\n0\n3\n{}\n".format(XY.shape[0]))
            for i in range(XY.shape[0]):
                fx, fy, fz = force_field[i]
                f.write(f"{i + 1} {fx} {fy} {fz}\n")
            f.write("$EndNodeData\n")
            
         # Marcar nós com apoio (restritos) com valor escalar para colorir
        if AP is not None and len(AP) > 0:
            apoio_marker = np.zeros(XY.shape[0])
            for nid, _, _ in AP:
                apoio_marker[int(nid) - 1] = 1.0  # Marca como restrito

            f.write("$NodeData\n1\n\"Apoios\"\n1\n0.0\n3\n0\n1\n{}\n".format(XY.shape[0]))
            for i, val in enumerate(apoio_marker):
                f.write(f"{i + 1} {val}\n")
            f.write("$EndNodeData\n")


def Exporta_para_Gmsh_GEO(filename, IJ, XY, U=None, TVM=None, detJ=None, AP=None, P=None, escala_setas=10.0):
    n_nodes = XY.shape[0]
    n_elements = IJ.shape[0]
    n_nos_por_elem = IJ.shape[1]
    
    if n_nos_por_elem == 3:
        gmsh_elem_type = 2
    elif n_nos_por_elem == 4:
        gmsh_elem_type = 3
    elif n_nos_por_elem == 8:
        gmsh_elem_type = 5
    else:
        raise ValueError("Número de nós por elemento não suportado.")

    node_ids = np.arange(1, n_nodes + 1)
    element_ids = np.arange(1, n_elements + 1)

    # Mapeamento das faces para hexaedro 8 nós
    face_map = [
        (0, 1, 2, 3),  # Face 1: bottom
        (4, 5, 6, 7),  # Face 2: top
        (0, 1, 5, 4),  # Face 3: front
        (1, 2, 6, 5),  # Face 4: right
        (2, 3, 7, 6),  # Face 5: back
        (3, 0, 4, 7)   # Face 6: left
    ]

    with open(filename, "w") as f:
        # Header
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        # Nós
        f.write(f"$Nodes\n{n_nodes}\n")
        for i, (x, y, z) in enumerate(XY, start=1):
            f.write(f"{i} {x} {y} {z}\n")
        f.write("$EndNodes\n")

        # Elementos
        f.write(f"$Elements\n{n_elements}\n")
        for eid, conn in zip(element_ids, IJ):
            conn_str = " ".join(str(nid) for nid in conn)
            f.write(f"{eid} {gmsh_elem_type} 2 0 1 {conn_str}\n")
        f.write("$EndElements\n")

        # Deslocamentos
        if U is not None:
            ux = U[0::3]
            uy = U[1::3]
            uz = U[2::3]

            f.write("$NodeData\n1\n\"Deslocamento Ux\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
            for nid, val in zip(node_ids, ux):
                f.write(f"{nid} {val}\n")
            f.write("$EndNodeData\n")

            f.write("$NodeData\n1\n\"Deslocamento Uy\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
            for nid, val in zip(node_ids, uy):
                f.write(f"{nid} {val}\n")
            f.write("$EndNodeData\n")

            f.write("$NodeData\n1\n\"Deslocamento Uz\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
            for nid, val in zip(node_ids, uz):
                f.write(f"{nid} {val}\n")
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
            for eid, val in zip(element_ids, detJ):
                f.write(f"{eid} {val}\n")
            f.write("$EndElementData\n")

        # Apoios (nós restritos)
        if AP is not None and len(AP) > 0:
            apoio_marker = np.zeros(n_nodes)
            for nid, _, _ in AP:
                apoio_marker[int(nid) - 1] = 1.0  # Marca nó restrito

            f.write("$NodeData\n1\n\"Apoios\"\n1\n0.0\n3\n0\n1\n{}\n".format(n_nodes))
            for i, val in enumerate(apoio_marker):
                f.write(f"{i + 1} {val}\n")
            f.write("$EndNodeData\n")

        # Forças aplicadas (P) distribuídas nos nós das faces (hexaedro)
        if P is not None and n_nos_por_elem == 8:
            force_field = np.zeros((n_nodes, 3))

            for eid, face_id, direcao, valor in P:
                conn = IJ[int(eid) - 1] - 1  # base zero
                face_nodes = [conn[i] for i in face_map[int(face_id) - 1]]
                for nid in face_nodes:
                    force_field[nid, direcao - 1] += (valor * escala_setas) / len(face_nodes)

            f.write("$NodeData\n1\n\"Carga Aplicada\"\n1\n0.0\n3\n0\n3\n{}\n".format(n_nodes))
            for i in range(n_nodes):
                fx, fy, fz = force_field[i]
                f.write(f"{i + 1} {fx} {fy} {fz}\n")
            f.write("$EndNodeData\n")

    # Agora gera o arquivo .geo para carregar o .msh e já configurar a visualização
    geo_filename = filename.replace(".msh", ".geo")
    with open(geo_filename, "w") as fgeo:
        fgeo.write(f'Merge "{filename}";\n\n')
        fgeo.write("// Configuração dos vetores\n")
        fgeo.write("View[3].VectorType = 5;\n")  # Flechas
        fgeo.write(f"View[3].MaxScale = {escala_setas};\n\n")

        fgeo.write("// Configuração do campo escalar Apoios\n")
        fgeo.write("View[4].ColorMap = 2;\n")  # Colormap quente (vermelho)
        fgeo.write("View[4].IntervalsType = 1;\n")
        fgeo.write("View[4].Min = 0;\n")
        fgeo.write("View[4].Max = 1;\n\n")

        fgeo.write("// Tamanho dos pontos\n")
        fgeo.write("General.PointSize = 8;\n")

    return geo_filename
    print(f"Arquivos '{filename}' e '{geo_filename}' gerados com sucesso.")


def AppendResultsToGmsh(filename_in, filename_out, U=None, TVM=None, detJ=None):
    with open(filename_in, "r") as f:
        lines = f.readlines()

    # Encontrar índice da linha do último "$EndElements"
    last_end_elements_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "$EndElements":
            last_end_elements_idx = i

    if last_end_elements_idx is None:
        raise ValueError("Arquivo msh inválido: não contém $EndElements")

    # Dividir arquivo em duas partes
    head = lines[:last_end_elements_idx + 1]
    tail = lines[last_end_elements_idx + 1:]

    # Extrair número de nós e elementos (assumindo formato padrão)
    n_nodes = 0
    n_elements = 0
    for i, line in enumerate(head):
        if line.strip() == "$Nodes":
            n_nodes = int(head[i+1].strip())
        if line.strip() == "$Elements":
            n_elements = int(head[i+1].strip())

    # Começa a montar as linhas com resultados
    results_lines = []

    if U is not None and n_nodes > 0:
        ux = U[0::3]
        uy = U[1::3]
        uz = U[2::3]

        def write_nodedata(name, values):
            results_lines.append("$NodeData\n")
            results_lines.append("1\n")  # número de conjuntos de dados
            results_lines.append(f"\"{name}\"\n")
            results_lines.append("1\n")  # time step ou algo assim
            results_lines.append("0.0\n")
            results_lines.append("3\n")  # dimension (3 = espaço)
            results_lines.append("0\n")
            results_lines.append("1\n")  # número de componentes por nó
            results_lines.append(f"{n_nodes}\n")
            for i, val in enumerate(values, 1):
                results_lines.append(f"{i} {val}\n")
            results_lines.append("$EndNodeData\n")

        write_nodedata("Deslocamento Ux", ux)
        write_nodedata("Deslocamento Uy", uy)
        write_nodedata("Deslocamento Uz", uz)

    if TVM is not None and n_elements > 0:
        results_lines.append("$ElementData\n")
        results_lines.append("1\n")
        results_lines.append("\"Tensão de von Mises\"\n")
        results_lines.append("1\n")
        results_lines.append("0.0\n")
        results_lines.append("3\n")
        results_lines.append("0\n")
        results_lines.append("1\n")
        results_lines.append(f"{n_elements}\n")
        for i, val in enumerate(TVM, 1):
            results_lines.append(f"{i} {val}\n")
        results_lines.append("$EndElementData\n")

    if detJ is not None and n_elements > 0:
        results_lines.append("$ElementData\n")
        results_lines.append("1\n")
        results_lines.append("\"Determinante do Jacobiano\"\n")
        results_lines.append("1\n")
        results_lines.append("0.0\n")
        results_lines.append("3\n")
        results_lines.append("0\n")
        results_lines.append("1\n")
        results_lines.append(f"{n_elements}\n")
        for i, val in enumerate(detJ, 1):
            results_lines.append(f"{i} {val}\n")
        results_lines.append("$EndElementData\n")

    # Escreve tudo no arquivo de saída
    with open(filename_out, "w") as f:
        f.writelines(head)
        f.writelines(results_lines)
        f.writelines(tail)

def Processar_malha3D(arquivo_msh, materiais, apoios, cargas):
    mesh = meshio.read(arquivo_msh)
     # 1. Coordenadas dos nós
    XY = mesh.points[:, :3]
    nn = XY.shape[0]

    # 2. Conectividade dos elementos
    if "hexahedron" in mesh.cells_dict:
        elements = mesh.cells_dict["hexahedron"]
        element_type = "hexahedron"
        IJ = elements + 1  # base 1
    elif "tetra" in mesh.cells_dict:
        elements = mesh.cells_dict["tetra"]
        element_type = "tetra"
        IJ = elements + 1
    else:
        raise ValueError("A malha deve conter elementos hexaédricos ou tetraédricos.")

    ne = IJ.shape[0]

    # 3. Materiais
    physical_tags = mesh.cell_data_dict["gmsh:physical"][element_type]
    tag_to_name = {v[0]: k for k, v in mesh.field_data.items()}
    MAT = []
    for tag in physical_tags:
        nome = tag_to_name.get(tag)
        if nome in materiais:
            E, nu = materiais[nome]
            MAT.append([E, nu])
        else:
            raise ValueError(f"Material '{nome}' não encontrado.")
    MAT = np.array(MAT)

    # 4. Apoios
    AP = []
    for name, tipo in apoios.items():
        for entity in mesh.field_data:
            if entity == name:
                tag, dim = mesh.field_data[entity]
                if dim == 2:
                    surfaces = mesh.cells_dict.get("quad", [])
                    phys_tags = mesh.cell_data_dict["gmsh:physical"]["quad"]
                    for i, t in enumerate(phys_tags):
                        if t == tag:
                            for node in surfaces[i]:
                                if "x"  in tipo:
                                    AP.append([node + 1, 1, 0])
                                if "y"  in tipo:
                                    AP.append([node + 1, 2, 0])
                                if "z"  in tipo:
                                    AP.append([node + 1, 3, 0])
                elif dim == 0:
                    for i, point in enumerate(mesh.cells_dict["vertex"]):
                        if mesh.cell_data_dict["gmsh:physical"]["vertex"][i] == tag:
                            if tipo in ("x", "engaste"):
                                AP.append([point[0] + 1, 1, 0])
                            if tipo in ("y", "engaste"):
                                AP.append([point[0] + 1, 2, 0])
                            if tipo in ("z", "engaste"):
                                    AP.append([node + 1, 3, 0])
    AP = np.unique(AP, axis=0)
    na = len(AP)

    # 5. Cargas
    P = []
    for nome, direcao, valor in cargas:
        tag, dim = mesh.field_data[nome]

        if dim == 2 and "quad" in mesh.cells_dict:
            quads = mesh.cells_dict["quad"]
            phys_tags = mesh.cell_data_dict["gmsh:physical"]["quad"]

            for i, t in enumerate(phys_tags):
                if t == tag:
                    face_nodes = set(quads[i])

                    for e, conn in enumerate(IJ - 1):  # conn = nós do hexaedro
                        for f, local_face in enumerate([
                            (0, 1, 2, 3),  # Face 1: Z = 0
                            (4, 5, 6, 7),  # Face 2: Z = t
                            (0, 1, 5, 4),  # Face 3: Y = 0
                            (1, 2, 6, 5),  # Face 4: X = L
                            (2, 3, 7, 6),  # Face 5: Y = h
                            (3, 0, 4, 7)   # Face 6: X = 0
                        ], 1):  # começa com 1

                            face_conn = {conn[i] for i in local_face}

                            if face_nodes == face_conn:
                                P.append([e + 1, f, direcao, valor])
                                break
        elif dim == 0:
            for i, point in enumerate(mesh.cells_dict["vertex"]):
                if mesh.cell_data_dict["gmsh:physical"]["vertex"][i] == tag:
                    if direcao == 1:
                        P.append([point[0] + 1, 1, 0, valor])
                    if direcao == 2:
                        P.append([point[0] + 1, 2, 0, valor])
                    if direcao == 3:
                        P.append([point[0] + 1, 3, 0, valor])
    P = np.array(P)
    nc = len(P)

    return nn, XY, ne, IJ, MAT, na, AP, nc, P


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
