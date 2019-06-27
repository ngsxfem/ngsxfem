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
    ngmesh = square.GenerateMesh(maxh=0.9, quad_dominated=False)
    mesh = Mesh (ngmesh)
    for i in range(lvl):
        mesh.Refine()
    return mesh


#trig version:
ngmesh = square.GenerateMesh(maxh=0.9, quad_dominated=False)
mesh = Mesh (ngmesh)
p1prolong = P1Prolongation(mesh)

prolongNeg = P1Prolongation(mesh)
prolongPos = P1Prolongation(mesh)



order = 1
nref = 2

Vh = H1(mesh, order = order, dirichlet=[], dgjumps = True)
lsetp1 = GridFunction(Vh)
InterpolateToP1(levelset,lsetp1)

def GetActiveDof(mesh,HASPOSNEG):
    Vh = H1(mesh, order = order, dirichlet=[], dgjumps = True)
    lsetp1 = GridFunction(Vh)
    InterpolateToP1(levelset,lsetp1)
    ci = CutInfo(mesh,lsetp1)
    return GetDofsOfElements(Vh,ci.GetElementsOfType(HASPOSNEG))

VhNeg = Compress( Vh, active_dofs = GetActiveDof(mesh, HASNEG) )
VhPos = Compress( Vh, active_dofs = GetActiveDof(mesh, HASPOS) )

CutVh = FESpace( [VhNeg, VhPos] )

CutProl = CompoundProlongation( CutVh )

CutProl.AddProlongation(prolongNeg)
CutProl.AddProlongation(prolongPos)

Vhc = Compress(Vh,active_dofs=GetActiveDof(mesh,HASNEG))
#p1prolong.Update(Vhc)
prolongNeg.Update(VhNeg)
prolongPos.Update(VhPos)
gfu0 = GridFunction( CutVh )

gfu0.components[0].Set( CoefficientFunction(x*x*y) )
gfu0.components[1].Set( CoefficientFunction(x*x*y+0.1) )
tmpvec = gfu0.vec.CreateVector()
for i in range(len(tmpvec)):
    tmpvec[i] = gfu0.vec[i]

cutgfu = GridFunction( CutVh )
cutgfu.vec[:] = 0.0

gfu = GridFunction(Vhc)
print("ndof",Vhc.ndof)
print( "dim vneg",CutVh.components[0].ndof )
print( "dim vpos",CutVh.components[1].ndof )
dofArr = [ 0 for i in range(nref+1) ]
dofArr[0] = CutVh.ndof
gf = []
gf.append(GridFunction(Compress(H1(mesh,order=1),active_dofs=GetActiveDof(mesh,HASNEG))))
for i in range(nref):
    mesh.Refine()
    Vh.Update()
    lsetp1.Update()
    InterpolateToP1(levelset,lsetp1)
    ci = CutInfo(mesh,lsetp1)
    Vhc.SetActiveDofs(GetActiveDof(mesh,HASNEG))
    Vhc.Update()

    CutVh.components[0].SetActiveDofs( GetActiveDof(mesh,HASNEG) )
    CutVh.components[1].SetActiveDofs( GetActiveDof(mesh,HASPOS) )
    CutVh.Update()
    CutProl.Update(CutVh)

    dofArr[i+1] = CutVh.ndof

    print( "dim vneg",CutVh.components[0].ndof )
    print( "dim vpos",CutVh.components[1].ndof )
    gfu.Update()
    cutgfu.Update()
    #p1prolong.Update(Vhc)

    
    nmesh = meshonlvl(i)
    gf.append(GridFunction(Compress(H1(nmesh,order=1),active_dofs=GetActiveDof(nmesh,HASNEG))))

gfdraw = GridFunction( CutVh )    

v = gfdraw.vec.CreateVector()    
w = gfdraw.vec.CreateVector()

Pv = gfdraw.vec.CreateVector()    
Rw = gfdraw.vec.CreateVector()

for i in range(len(v)):
    v[i] = i
for i in range(len(w)):
    w[i] = 1.0/(i+1)

Rw.data = w    
CutProl.Restrict(2,Rw)
CutProl.Restrict(1,Rw)

#v[dofArr[0]: ] = 0.0
Pv.data = v
CutProl.Prolongate(1,Pv)
CutProl.Prolongate(2,Pv)

# print(v)
# print(Pv)
# print(w)
# print(Rw)

print(InnerProduct(Pv,w))
print(InnerProduct(v,Rw))

input('w')

    
v = gfu.vec.CreateVector()    
w = gfu.vec.CreateVector()


prolvec = cutgfu.vec.CreateVector()
for i in range( len(tmpvec) ):
    prolvec[i] = tmpvec[i]


print('before prol')
for i in range(nref):
    CutProl.Prolongate(i+1,prolvec)


offset = [ 0, len(gfdraw.components[0].vec) ]
lens   = [ len(gfdraw.components[i].vec) for i in [0,1] ]

for i in [0,1]:
    for k in range(offset[i], offset[i]+lens[i] ):
        gfdraw.components[i].vec[k-offset[i]] = prolvec[k]

u_coef = IfPos(lsetp1, gfdraw.components[1], gfdraw.components[0])
Draw(u_coef,mesh,"u")


input('w')

for i in range(len(v)):
    v[i] = i
w.data = v
# Pv = gfu.vec.CreateVector()    
# Rw = gfu.vec.CreateVector()


for i in range(len(gf[0].vec)):
    gf[0].vec[i] = i


for i in range(nref+1):
    if i > 0:
       p1prolong.Prolongate(i,w) 
    #print(w)
    nmesh = meshonlvl(i)
    gfvis = GridFunction(Compress(H1(nmesh,order=1),active_dofs=GetActiveDof(nmesh,HASNEG)))
    for i in range(len(gfvis.vec)):
        gfvis.vec[i] = w[i]
    Draw(gfvis,nmesh,"gf")
    input(i)
    
    
# v = gfu.vec.CreateVector()    
# w = gfu.vec.CreateVector()

# Pv = gfu.vec.CreateVector()    
# Rw = gfu.vec.CreateVector()

# for i in range(len(v)):
#     v[i] = i
# for i in range(len(w)):
#     w[i] = 1.0/(i+1)

# Rw.data = w    
# p1prolong.Restrict(2,Rw)
# p1prolong.Restrict(1,Rw)

# Pv.data = v
# p1prolong.Prolongate(1,Pv)
# p1prolong.Prolongate(2,Pv)

# print(InnerProduct(Pv,w))
# print(InnerProduct(v,Rw))
    
# # gfu.vec[:] = 1.0 

# # input("")
# print(gfu.vec)

# p1prolong.Restrict(2,gfu.vec)
# print(gfu.vec)
# input("")
# p1prolong.Restrict(1,gfu.vec)

# print(gfu.vec)
# input("")

