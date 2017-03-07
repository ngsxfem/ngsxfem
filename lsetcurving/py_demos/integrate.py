# integration on lset domains
from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

#square domain [-1,1]x[-1,1]
def Make2DProblem(maxh=2):
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([-1,-1],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=False))
    return mesh;

# circle with radius 0.5
levelset = sqrt(x*x+y*y)-0.5
referencevals = { POS : 4.0-0.25*pi, NEG : 0.25*pi, IF : pi }

mesh = Make2DProblem(maxh=0.5)

order = 5
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)

errors_uncurved = dict()
errors_curved = dict()
eoc_uncurved = dict()
eoc_curved = dict()

for key in [NEG,POS,IF]:
    errors_curved[key] = []
    errors_uncurved[key] = []
    eoc_curved[key] = []
    eoc_uncurved[key] = []


for reflevel in range(6):

    if(reflevel > 0):
        mesh.Refine()
    
    f = CoefficientFunction (1.0)

    for key in [NEG,POS,IF]:
        integrals_uncurved = Integrate(levelset_domain = { "levelset" : levelset, "domain_type" : key},
                                       cf=f, mesh=mesh, order = order)
        # Applying the mesh deformation
        deformation = lsetmeshadap.CalcDeformation(levelset)
        mesh.SetDeformation(deformation)
        integrals_curved = Integrate(levelset_domain = { "levelset" : levelset, "domain_type" : key},
                                       cf=f, mesh=mesh, order = order)
        # Unapply the mesh deformation (for refinement)
        mesh.UnsetDeformation()

        errors_curved[key].append(abs(integrals_curved - referencevals[key]))
        errors_uncurved[key].append(abs(integrals_uncurved - referencevals[key]))
    # refine cut elements:
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)

for key in [NEG,POS,IF]:
    eoc_curved[key] = [log(a/b)/log(2) for (a,b) in zip (errors_curved[key][0:-1],errors_curved[key][1:]) ]
    eoc_uncurved[key] = [log(a/b)/log(2) for (a,b) in zip (errors_uncurved[key][0:-1],errors_uncurved[key][1:]) ]

print("errors (uncurved):  \n{}\n".format(errors_uncurved))
print("   eoc (uncurved):  \n{}\n".format(   eoc_uncurved))
print("errors (  curved):  \n{}\n".format(  errors_curved))
print("   eoc (  curved):  \n{}\n".format(     eoc_curved))

Draw(levelset,mesh,"levelset")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")


