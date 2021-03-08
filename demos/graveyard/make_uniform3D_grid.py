"""
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.meshing import Element2D, Element3D, MeshPoint, FaceDescriptor
from netgen.meshing import Mesh as NetMesh
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from ngsolve import Mesh as NGSMesh


# ----------------------------------- MAIN ------------------------------------
def MakeUniform3DGrid(quads=False, N=5, P1=(0, 0, 0), P2=(1, 1, 1)):

    if not quads:
        raise BaseException("Only hex cube so far...")

    Lx = P2[0] - P1[0]
    Ly = P2[1] - P1[1]
    Lz = P2[2] - P1[2]

    cube = OrthoBrick(Pnt(P1[0], P1[1], P1[2]), Pnt(P2[0], P2[1], P2[2])).bc(1)
    geom = CSGeometry()
    geom.Add(cube)
    netmesh = NetMesh()
    netmesh.SetGeometry(geom)
    netmesh.dim = 3

    pids = []
    for i in range(N + 1):
        for j in range(N + 1):
            for k in range(N + 1):
                pids.append(netmesh.Add(MeshPoint(Pnt(P1[0] + Lx * i / N,
                                                      P1[1] + Ly * j / N,
                                                      P1[2] + Lz * k / N))))

    for i in range(N):
        for j in range(N):
            for k in range(N):
                base = i * (N + 1) * (N + 1) + j * (N + 1) + k
                baseup = base + (N + 1) * (N + 1)
                pnum = [base, base + 1, base + (N + 1) + 1, base + (N + 1),
                        baseup, baseup + 1, baseup + (N + 1) + 1,
                        baseup + (N + 1)]
                elpids = [pids[p] for p in pnum]
                netmesh.Add(Element3D(1, elpids))

    def AddSurfEls(p1, dx, dy, facenr):
        for i in range(N):
            for j in range(N):
                base = p1 + i * dx + j * dy
                pnum = [base, base + dx, base + dx + dy, base + dy]
                elpids = [pids[p] for p in pnum]
                netmesh.Add(Element2D(facenr, elpids))

    netmesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    AddSurfEls(0, 1, N + 1, facenr=1)

    netmesh.Add(FaceDescriptor(surfnr=2, domin=1, bc=1))
    AddSurfEls(0, (N + 1) * (N + 1), 1, facenr=2)

    netmesh.Add(FaceDescriptor(surfnr=3, domin=1, bc=1))
    AddSurfEls(0, N + 1, (N + 1) * (N + 1), facenr=3)

    netmesh.Add(FaceDescriptor(surfnr=4, domin=1, bc=1))
    AddSurfEls((N + 1)**3 - 1, -(N + 1), -1, facenr=1)

    netmesh.Add(FaceDescriptor(surfnr=5, domin=1, bc=1))
    AddSurfEls((N + 1)**3 - 1, -(N + 1) * (N + 1), -(N + 1), facenr=1)

    netmesh.Add(FaceDescriptor(surfnr=6, domin=1, bc=1))
    AddSurfEls((N + 1)**3 - 1, -1, -(N + 1) * (N + 1), facenr=1)

    netmesh.Compress()
    mesh = NGSMesh(netmesh)
    mesh.ngmesh.Save("tmp.vol.gz")
    mesh = NGSMesh("tmp.vol.gz")
    return mesh
