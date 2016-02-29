# solve the Laplace Beltrami equation on a sphere
# with Dirichlet boundary condition u = 0
from math import pi

# ngsolve stuff
from ngsolve import *

# xfem and trace fem stuff
import libngsxfem_py.xfem as xfem                                 

# geometry stuff
# unit square [-1.41,1.41]^2 (boundary index 1 everywhere)
def MakeSquare():
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([-1.41,-1.41],[1.41,1.41],bc=1)
    return square;

# unit cube [-1.41,1.41]^3 (boundary index 1 everywhere)
def MakeCube():
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    cube = CSGeometry()
    cube.Add (OrthoBrick(Pnt(-1.41,-1.41,-1.41), Pnt(1.41,1.41,1.41)))
    return cube;

# meshing
from netgen.meshing import MeshingParameters
ngsglobals.msg_level = 1
mesh = Mesh (MakeCube().GenerateMesh(maxh=4.0, quad_dominated=False))
mesh.Refine()
mesh.Refine()

# The mesh deformation calculator
from xfem.lsetcurv import *
order = 2
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.5, discontinuous_qn=True)

# The level set
levelset = sqrt(x*x+y*y+z*z) - 1

# the solution
sol = sin(pi*z)

from xfem.tracefem import *
Vh = H1(mesh, order=order, dirichlet=[])
Vh_tr = TraceFESpace(mesh, Vh, levelset)

a = BilinearForm(Vh_tr, symmetric = True, flags = { })
a += TraceMass(1.0)
a += TraceLaplaceBeltrami(1.0)

f = LinearForm(Vh_tr)
f += TraceSource(sin(pi*z)*(pi*pi*(1-z*z)+1)+cos(pi*z)*2*pi*z)

c = Preconditioner(a, type="local", flags= { "test" : True })
#c = Preconditioner(a, type="direct", flags= { "inverse" : "pardiso", "test" : True })

u = GridFunction(Vh_tr)

vtk = VTKOutput(ma=mesh,coefs=[levelset,lsetmeshadap.lset_ho,lsetmeshadap.lset_p1,lsetmeshadap.deform,u,sol],names=["lset","lsetho","lsetp1","deform","u","uexact"],filename="vtkout_simple_",subdivision=2)

def SolveProblem():
    # Calculation of the deformation:
    deformation = lsetmeshadap.CalcDeformation(levelset)
    # Applying the mesh deformation
    mesh.SetDeformation(deformation)

    Vh.Update()
    Vh_tr.Update()
    u.Update()
    
    a.Assemble();
    f.Assemble();
    c.Update();
    
    solvea = CGSolver( mat=a.mat, pre=c.mat, complex=False, printrates=False, precision=1e-8, maxsteps=2000)
    # # the boundary value problem to be solved on each level
    u.vec.data = solvea * f.vec;
    
    print ("Number of iterations:", solvea.GetSteps())

def PostProcess():
    print ("Max. dist of discr. intf. to exact level set fct:", lsetmeshadap.CalcMaxDistance(levelset))
    print ("Max. dist of discr. intf. to discr. level set f.:", lsetmeshadap.CalcMaxDistance(lsetmeshadap.lset_ho))
    CalcTraceDiff( u, sol, intorder=2*order);
    mesh.UnsetDeformation()
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)
    vtk.Do()
    Draw(lsetmeshadap.lset_p1,mesh,"lsetp1")
    Draw(lsetmeshadap.deform,mesh,"deformation")
    Draw(u,mesh,"u",draw_surf=False)
    
SolveProblem()
PostProcess()
for i in range(2):
    mesh.Refine()
    SolveProblem()
    PostProcess()

    
