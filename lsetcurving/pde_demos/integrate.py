# solve the Laplace Beltrami equation on a sphere
# with Dirichlet boundary condition u = 0
from math import pi
# ngsolve stuff
from ngsolve import *
# xfem and trace fem stuff
import libngsxfem_py.xfem as xfem                                 
# basic xfem functionality
from xfem.basics import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *

# 2D: circle configuration
def Make2DProblem(maxh=2):
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([-1,-1],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=False))
    return mesh;

levelset = sqrt(x*x+y*y)-0.5
referencevals = [4.0-0.25*pi,0.25*pi]

mesh = Make2DProblem(maxh=0.5)

order = 5 # Pk+1 Pk
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)

errors_uncurved = []
errors_curved = []

for reflevel in range(6):

    if(reflevel > 0):
        mesh.Refine()
    
    testcoef = CoefficientFunction (1.0)
    integrals_uncurved = IntegrateX(testcoef,levelset,mesh,order=order)

    # Applying the mesh deformation
    deformation = lsetmeshadap.CalcDeformation(levelset)
    mesh.SetDeformation(deformation)

    integrals_curved = IntegrateX(testcoef,levelset,mesh,order=order)

    # Unapply the mesh deformation (for refinement)
    mesh.UnsetDeformation()

    # refine cut elements:
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)
    
    errors_uncurved.append([ abs(a-b) for (a,b) in zip (integrals_uncurved,referencevals)])
    errors_curved.append([ abs(a-b) for (a,b) in zip (integrals_curved,referencevals)])
    
eoc_uncurved = [[ log(a[j]/b[j])/log(2) for j in [0,1] ] for (a,b) in zip (errors_uncurved[0:-1],errors_uncurved[1:]) ]
eoc_curved = [[ log(a[j]/b[j])/log(2) for j in [0,1] ] for (a,b) in zip (errors_curved[0:-1],errors_curved[1:]) ]

print("errors (uncurved):  {}\n".format(errors_uncurved))
print("   eoc (uncurved):  {}\n".format(   eoc_uncurved))
print("errors (  curved):  {}\n".format(  errors_curved))
print("   eoc (  curved):  {}\n".format(     eoc_curved))

Draw(levelset,mesh,"levelset")
Draw(lsetmeshadap.deform,mesh,"deformation")
Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")


