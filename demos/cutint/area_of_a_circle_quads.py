from math import pi
from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
#mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))
r=0.6

#domains = [NEG,POS,IF]
domains = [NEG, POS, IF]

levelset = sqrt(x*x+y*y)-r
referencevals = { POS : 1-pi*r*r/4, NEG : pi*r*r/4, IF : r*pi/2}

n_ref = 8
order = 2
errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []


lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)

for i in range(n_ref):
  V = H1(mesh,order=1)
  lset_approx = GridFunction(V)
  InterpolateToP1(levelset,lset_approx)
  #Draw(lset_approx,mesh,"lset_p1")

  f = CoefficientFunction(1)

  deformation = lsetmeshadap.CalcDeformation(levelset)
  #mesh.SetDeformation(deformation)
  #Draw(deformation,mesh,"deformation")

  for key in domains:
    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : key},
                                     cf=f, mesh=mesh, order = order)
    print("Result of Integration Reflevel ",i,", Key ",key," : ", integral)
    errors[key].append(abs(integral - referencevals[key]))

  mesh.UnsetDeformation()

  if i < n_ref - 1:
    # RefineAtLevelSet(gf=lset_approx)
    mesh.Refine()

for key in domains:
  eoc[key] = [log(errors[key][i+1]/errors[key][i])/log(0.5) for i in range(n_ref-1)]

print("L2-errors:", errors)
#l2_eoc = [log(errors[i+1]/errors[i])/log(0.5) for i in range(n_ref-1)]
print("experimental order of convergence (L2):", eoc)

for key in domains:
    for i in range(2, len(eoc[key])):
        eoc_v = eoc[key][i]
        if eoc_v < 1.8:
            raise RuntimeError("Order of Convergence < 1.8 detected!")
