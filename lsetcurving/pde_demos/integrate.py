# integration on lset domains
from math import pi
# ngsolve stuff
from ngsolve import *
# xfem stuff
import libngsxfem_py.xfem as xfem                                 
# basic xfem functionality
from xfem.basics import *
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
referencevals = { "posdomain" : 4.0-0.25*pi, "negdomain" : 0.25*pi, "interface" : pi }

mesh = Make2DProblem(maxh=0.5)

order = 5
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)

domains = interface_and_volume_domains

errors_uncurved = dict()
errors_curved = dict()
eoc_uncurved = dict()
eoc_curved = dict()
for key in domains:
    errors_curved[key] = []
    errors_uncurved[key] = []
    eoc_curved[key] = []
    eoc_uncurved[key] = []


for reflevel in range(6):

    if(reflevel > 0):
        mesh.Refine()
    
    f = CoefficientFunction (1.0)
    integrals_uncurved = IntegrateX(levelset,mesh,cf_neg=f,cf_pos=f,cf_interface=f,order=order, domains=domains)
    # Applying the mesh deformation
    deformation = lsetmeshadap.CalcDeformation(levelset)
    mesh.SetDeformation(deformation)

    integrals_curved = IntegrateX(levelset,mesh,cf_neg=f,cf_pos=f,cf_interface=f,order=order, domains=domains)
    # Unapply the mesh deformation (for refinement)
    mesh.UnsetDeformation()

    # refine cut elements:
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)

    for key in domains:
        errors_curved[key].append(abs(integrals_curved[key] - referencevals[key]))
        errors_uncurved[key].append(abs(integrals_uncurved[key] - referencevals[key]))

for key in domains:
    eoc_curved[key] = [log(a/b)/log(2) for (a,b) in zip (errors_curved[key][0:-1],errors_curved[key][1:]) ]
    eoc_uncurved[key] = [log(a/b)/log(2) for (a,b) in zip (errors_uncurved[key][0:-1],errors_uncurved[key][1:]) ]

print("errors (uncurved):  \n{}\n".format(errors_uncurved))
print("   eoc (uncurved):  \n{}\n".format(   eoc_uncurved))
print("errors (  curved):  \n{}\n".format(  errors_curved))
print("   eoc (  curved):  \n{}\n".format(     eoc_curved))

Draw(levelset,mesh,"levelset")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")


