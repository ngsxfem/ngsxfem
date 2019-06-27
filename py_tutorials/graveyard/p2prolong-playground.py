from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters
from ngsolve import *
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

r2 = 3/4 # outer radius
r1 = 1/4 # inner radius
rc = (r1+r2) / 2.0
rr = (r2-r1) / 2.0
r = sqrt(x*x+y*y)
levelset = IfPos(r-rc, r-rc-rr,rc-r-rr)

# Geometry 
square = SplineGeometry()
square.AddRectangle([-1,-1],[1,1],bc=1)

def meshonlvl(lvl):
    ngmesh = square.GenerateMesh(maxh=0.1, quad_dominated=False)
    mesh = Mesh (ngmesh)
    for i in range(lvl):
        mesh.Refine()
    return mesh


#trig version:
ngmesh = square.GenerateMesh(maxh=0.1, quad_dominated=False)
mesh = Mesh (ngmesh)
p1prolong = P2CutProlongation(mesh)

polorder = 2

order = 1
nref = 2

Vh = H1(mesh, order = polorder, dirichlet=[], dgjumps = True)
lsetp1 = GridFunction(Vh)
InterpolateToP1(levelset,lsetp1)

def GetActiveDof(mesh):
    Vh = H1(mesh, order = polorder, dirichlet=[], dgjumps = True)
    lsetp1 = GridFunction(Vh)
    InterpolateToP1(levelset,lsetp1)
    ci = CutInfo(mesh,lsetp1)
    return GetDofsOfElements(Vh,ci.GetElementsOfType(HASNEG))

Vhc = Compress(Vh,active_dofs=GetActiveDof(mesh))
p1prolong.Update(Vhc)
gfu = GridFunction(Vhc)
print("ndof",Vhc.ndof)
gf = []
gf.append(GridFunction(Compress(H1(mesh,order=polorder),active_dofs=GetActiveDof(mesh))))
gf[0].Set(x*x*x*y)
vc = gf[0].vec.CreateVector()
vc.data = gf[0].vec

for i in range(nref):
    mesh.Refine()
    Vh.Update()
    lsetp1.Update()
    InterpolateToP1(levelset,lsetp1)
    ci = CutInfo(mesh,lsetp1)
    Vhc.SetActiveDofs(GetActiveDof(mesh))
    Vhc.Update()
    print(Vh.ndof)
    print(Vhc.ndof)
    gfu.Update()
    p1prolong.Update(Vhc)

    
    nmesh = meshonlvl(i)
    gf.append(GridFunction(Compress(H1(nmesh,order=polorder),active_dofs=GetActiveDof(nmesh))))
    
input('mesh updated')
    
v = gfu.vec.CreateVector()    
w = gfu.vec.CreateVector()

for i in range(len(vc)):
    v[i] = vc[i]
w.data = v
# Pv = gfu.vec.CreateVector()    
# Rw = gfu.vec.CreateVector()


# for i in range(len(gf[0].vec)):
#     gf[0].vec[i] = i


for i in range(nref+1):
    if i > 0:
       p1prolong.Prolongate(i,w) 
    #print(w)
    nmesh = meshonlvl(i)
    gfvis = GridFunction(Compress(H1(nmesh,order=polorder),active_dofs=GetActiveDof(nmesh)))
    for i in range(len(gfvis.vec)):
        gfvis.vec[i] = w[i]
    Draw(gfvis,nmesh,"gf")
    input(i)
    
    
v = gfu.vec.CreateVector()    
w = gfu.vec.CreateVector()

Pv = gfu.vec.CreateVector()    
Rw = gfu.vec.CreateVector()

for i in range(len(v)):
    v[i] = i
for i in range(len(w)):
    w[i] = 1.0/(i+1)

Rw.data = w    
p1prolong.Restrict(2,Rw)
p1prolong.Restrict(1,Rw)

Pv.data = v
p1prolong.Prolongate(1,Pv)
p1prolong.Prolongate(2,Pv)

print(InnerProduct(Pv,w))
print(InnerProduct(v,Rw))
    
# # gfu.vec[:] = 1.0 

# # input("")
# print(gfu.vec)

# p1prolong.Restrict(2,gfu.vec)
# print(gfu.vec)
# input("")
# p1prolong.Restrict(1,gfu.vec)

# print(gfu.vec)
# input("")

