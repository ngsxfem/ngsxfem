from netgen.geom2d import unit_square, MakeCircle, SplineGeometry, MakeRectangle
from netgen.meshing import Element0D, Element1D, Element2D, MeshPoint, FaceDescriptor, Mesh as NetMesh
from netgen.csg import Pnt
from ngsolve import Mesh as NGSMesh

def MakeUniform2DGrid(quads = False, N=5, P1=(0,0),P2=(1,1)):
  geom = SplineGeometry()
  geom.AddRectangle(P1,P2,bc=1)

  Lx = P2[0]-P1[0]
  Ly = P2[1]-P1[1]
  
  netmesh = NetMesh()
  netmesh.SetGeometry(geom)
  netmesh.dim = 2
  pnums = []
  for i in range(N + 1):
      for j in range(N + 1):
          pnums.append(netmesh.Add(MeshPoint(Pnt(P1[0] + Lx * i / N, P1[1] + Ly * j / N, 0))))

  netmesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
  netmesh.SetMaterial(1, "mat")
  for j in range(N):
      for i in range(N):
          if quads:
              netmesh.Add(Element2D(1, [pnums[i + j * (N + 1)], pnums[i + (j + 1) * (N + 1)], pnums[i + 1 + (j + 1) * (N + 1)], pnums[i + 1 + j * (N + 1)]]))
          else:
              netmesh.Add(Element2D(1, [pnums[i + j * (N + 1)], pnums[i + (j + 1) * (N + 1)], pnums[i + 1 + j * (N + 1)]]))
              netmesh.Add(Element2D(1, [pnums[i + (j + 1) * (N + 1)], pnums[i + 1 + (j + 1) * (N + 1)], pnums[i + 1 + j * (N + 1)]]))
              
  
  for j in range(N):
      netmesh.Add(Element1D([pnums[N + j * (N + 1)], pnums[N + (j + 1) * (N + 1)]], index=1))
      netmesh.Add(Element1D([pnums[0 + j * (N + 1)], pnums[0 + (j + 1) * (N + 1)]], index=1))
  
      
  for i in range(N):
      netmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=1))
      netmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=1))
  
  netmesh.SetGeometry(geom)
  netmesh.Compress()
  mesh = NGSMesh(netmesh)
  mesh.ngmesh.Save("tmp.vol.gz")
  mesh = NGSMesh("tmp.vol.gz")
  return mesh 
