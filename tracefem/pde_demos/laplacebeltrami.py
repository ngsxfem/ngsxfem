# solve the Laplace Beltrami equation on a sphere
# with Dirichlet boundary condition u = 0
from math import pi

# ngsolve stuff
from ngsolve import *

# xfem and trace fem stuff
import libngsxfem_py.xfem as xfem                                 

# for plotting (convergence plots)
import matplotlib.pyplot as plt

# for asking for interactive shell or not
import sys    


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
# criss cross mesh with midpoint
mesh.Refine()

# The mesh deformation calculator
from xfem.lsetcurv import *
order = 2
lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)

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

#vtk = VTKOutput(ma=mesh,coefs=[levelset,lsetmeshadap.lset_ho,lsetmeshadap.lset_p1,lsetmeshadap.deform,u,sol],names=["lset","lsetho","lsetp1","deform","u","uexact"],filename="vtkout_simple_",subdivision=2)

global last_num_its
last_num_its = 0
statistics = []
level = 0

def SolveProblem():
    # Calculation of the deformation:
    deformation = lsetmeshadap.CalcDeformation(levelset)
    # Applying the mesh deformation
    mesh.SetDeformation(deformation)

    Vh.Update()
    # print("StdFESpace NDof:", Vh.ndof)
    Vh_tr.Update()
    u.Update()
    
    a.Assemble();
    f.Assemble();
    c.Update();
    
    solvea = CGSolver( mat=a.mat, pre=c.mat, complex=False, printrates=False, precision=1e-8, maxsteps=2000)
    # # the boundary value problem to be solved on each level
    u.vec.data = solvea * f.vec;
    
    global last_num_its
    last_num_its = solvea.GetSteps()
    mesh.UnsetDeformation()

def PostProcess():
    maxdistlset = lsetmeshadap.CalcMaxDistance(levelset);
    # maxdistlsetho = lsetmeshadap.CalcMaxDistance(lsetmeshadap.lset_ho);
    mesh.SetDeformation(lsetmeshadap.deform)
    [l2diff,maxdiff] = CalcTraceDiff( u, sol, intorder=2*order+2)
    mesh.UnsetDeformation()

    print ("The mesh Refinement level :", level)
    print ("Number of d.o.f. Std. -FES:", Vh.ndof)
    print ("Number of d.o.f. Trace-FES:", Vh_tr.ndof)
    print ("Max. distance to interface:", maxdistlset)
  # print ("Max. distance to ho-intfce:", maxdistlsetho)
    print ("L2-norm error on interface:", l2diff)
    print ("maxnorm error on interface:", maxdiff)
    print ("Numb. of CG-Solver iterat.:", last_num_its)
    statistics.append ( (level, Vh_tr.ndof, maxdistlset, l2diff, maxdiff, last_num_its) )
    
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)
    #vtk.Do()
    Draw(lsetmeshadap.lset_p1,mesh,"lsetp1")
    Draw(lsetmeshadap.deform,mesh,"deformation")
    Draw(u,mesh,"u",draw_surf=False)


def PlotConvergence():
    lvl,ndof,geomerr,l2err,maxerr,itcnts = zip(*statistics)


    plt.figure(0)
    plt.yscale('log')
    plt.xlabel("level")
    
    plt.plot(lvl,l2err, "-*")
    plt.plot(geomerr, "-+")
    plt.legend(["L2error","geometry (max) error"])

    # plt.figure(1)
    # plt.yscale('log')
    # plt.xlabel("level")

    # plt.plot(lvl,ndof, "-*")
    # plt.plot(itcnts, "-+")
    # plt.legend(["D.o.f.","iterations"])

    # plt.figure(2)
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel("ndof")
    # plt.ylabel("error")

    # plt.plot(ndof,l2err, "-*")
    # plt.legend(["L2error"])
    
    plt.ion()
    plt.show()

    
    l2conv = [ log(l2err[i]/l2err[i-1])/log(0.5) for i in range(1,len(l2err))]
    maxconv = [ log(maxerr[i]/maxerr[i-1])/log(0.5) for i in range(1,len(maxerr))]
    geomconv = [ log(geomerr[i]/geomerr[i-1])/log(0.5) for i in range(1,len(geomerr))]
    print ("l2err convergence orders (eoc):", l2conv)
    print ("geom. convergence orders (eoc):", geomconv)
    if (not hasattr(sys,'ps1')):
        input("<press enter to quit>")

def InitialSolve():
    with TaskManager():
        SolveProblem()
    PostProcess()
    
def RefineAndSolve():
    with TaskManager():
        mesh.Refine()
    global level
    level = level + 1
    with TaskManager():
        SolveProblem()
    PostProcess()

def LaplaceBeltramiOnSphere(reflvls = 1):
    InitialSolve()
    for i in range(reflvls-1):
        RefineAndSolve()
    if (reflvls > 1):
        PlotConvergence()
    
if __name__ == "__main__":
    LaplaceBeltramiOnSphere(reflvls = 5)
