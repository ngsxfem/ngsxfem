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
# for making a directory (if it doesn't exist)
import os
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *
# For TraceFEM-Integrators (convenience)
from xfem.tracefem import *

# 2D: circle configuration
def Make2DProblem():
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([-1.41,-1.41],[1.41,1.41],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=1, quad_dominated=False))
    mesh.Refine()

    problem = {"Diffusion" : 1.0,
               "Convection" : None,
               "Reaction" : 1.0,
               "Source" : sin(pi*y)*(1+pi*pi*(1-y*y))+pi*y*cos(pi*y),
               "Solution" : sin(pi*y),
               "Levelset" : sqrt(x*x+y*y)-1,
               "Mesh" : mesh
              }
    return problem;

# 3D: circle configuration
def Make3DProblem():
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    cube = CSGeometry()
    cube.Add (OrthoBrick(Pnt(-1.41,-1.41,-1.41), Pnt(1.41,1.41,1.41)))
    mesh = Mesh (cube.GenerateMesh(maxh=1, quad_dominated=False))
    mesh.Refine()

    problem = {"Diffusion" : 1.0,
               "Convection" : None,
               "Reaction" : 1.0,
               "Source" : sin(pi*z)*(pi*pi*(1-z*z)+1)+cos(pi*z)*2*pi*z,
               "Solution" : sin(pi*z),
               "Levelset" : sqrt(x*x+y*y+z*z)-1,
               "Mesh" : mesh
              }
    return problem;

problemdata = Make3DProblem()
mesh = problemdata["Mesh"]
order = 2

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)

### Setting up discrete variational problem

Vh = H1(mesh, order=order, dirichlet=[])
Vh_tr = TraceFESpace(mesh, Vh, problemdata["Levelset"])

a = BilinearForm(Vh_tr, symmetric = True, flags = {"eliminate_internal" : False})
if (problemdata["Reaction"] != None):
    a += TraceMass(problemdata["Reaction"])
if (problemdata["Diffusion"] != None):
    a += TraceLaplaceBeltrami(problemdata["Diffusion"])
    a += NormalLaplaceStabilization(problemdata["Diffusion"],lsetmeshadap.lset_p1.Deriv())
                                    #/(sqrt(lsetmeshadap.lset_p1.Deriv()*lsetmeshadap.lset_p1.Deriv())))
if (problemdata["Convection"] != None):
    a += TraceConvection(problemdata["Convection"])

f = LinearForm(Vh_tr)
if (problemdata["Source"] != None):
    f += TraceSource(problemdata["Source"])

c = Preconditioner(a, type="local", flags= { "test" : True })
#c = Preconditioner(a, type="direct", flags= { "inverse" : "pardiso", "test" : True })

u = GridFunction(Vh_tr)


if not os.path.exists("LapBeltrResults"):
    os.makedirs("LapBeltrResults")
vtk = VTKOutput(ma=mesh,coefs=[problemdata["Levelset"],lsetmeshadap.lset_ho,lsetmeshadap.lset_p1,lsetmeshadap.deform,u,problemdata["Solution"]],names=["lset","lsetho","lsetp1","deform","u","uexact"],filename="LapBeltrResults/vtkout_",subdivision=2)

global last_num_its
last_num_its = 0
statistics = []
level = 0

def SolveProblem():
    # Calculation of the deformation:
    deformation = lsetmeshadap.CalcDeformation(problemdata["Levelset"])
    # Applying the mesh deformation
    mesh.SetDeformation(deformation)

    Vh.Update()
    # print("StdFESpace NDof:", Vh.ndof)
    Vh_tr.Update()
    u.Update()
    
    a.Assemble();
    f.Assemble();
    c.Update();
    
    solvea = CGSolver( mat=a.mat, pre=c.mat, complex=False, printrates=False, precision=1e-8, maxsteps=200000)
    # # the boundary value problem to be solved on each level
    u.vec.data = solvea * f.vec;
    
    global last_num_its
    last_num_its = solvea.GetSteps()
    mesh.UnsetDeformation()

def PostProcess(vtkout = False):
    maxdistlset = lsetmeshadap.CalcMaxDistance(problemdata["Levelset"]);
    # maxdistlsetho = lsetmeshadap.CalcMaxDistance(lsetmeshadap.lset_ho);
    mesh.SetDeformation(lsetmeshadap.deform)
    [l2diff,maxdiff] = CalcTraceDiff( u, problemdata["Solution"], intorder=2*order+2)
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
    if (vtkout):
        vtk.Do()
    Draw(lsetmeshadap.lset_p1,mesh,"lsetp1")
    Draw(lsetmeshadap.deform,mesh,"deformation")
    Draw(u,mesh,"u",draw_surf=False)


def ConvergenceStudy(plot=True):
    lvl,ndof,geomerr,l2err,maxerr,itcnts = zip(*statistics)

    fo = open("LapBeltrResults/l2err.csv","w"); fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(lvl,l2err)])
    fo = open("LapBeltrResults/ndof.csv","w"); fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(lvl,ndof)])
    fo = open("LapBeltrResults/geomerr.csv","w"); fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(lvl,geomerr)])
    fo = open("LapBeltrResults/maxerr.csv","w"); fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(lvl,maxerr)])
    fo = open("LapBeltrResults/itcnts.csv","w"); fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(lvl,itcnts)])

    if (plot):
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
    if (plot and not hasattr(sys,'ps1')):
        input("<press enter to quit>")

def InitialSolve(vtkout = False):
    with TaskManager():
        SolveProblem()
    PostProcess(vtkout = vtkout)
    
def RefineAndSolve(vtkout = False):
    with TaskManager():
        mesh.Refine()
    global level
    level = level + 1
    with TaskManager():
        SolveProblem()
    PostProcess(vtkout = vtkout)

def LaplaceBeltramiOnSphere(reflvls = 1, vtkout = False, plot = True):
    InitialSolve(vtkout = vtkout)
    for i in range(reflvls-1):
        RefineAndSolve(vtkout = vtkout)
    if (reflvls > 1):
        ConvergenceStudy(plot = plot)


import argparse        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Solve Laplace-Betrami example problem.')
    parser.add_argument('--reflvls', metavar='N', type=int, default=2, help='number of refinement levels (>=1)')
    parser.add_argument('--vtkout', dest='vtkout', action='store_true', help='enable VTKOutput')
    parser.add_argument('--no-vtkout', dest='vtkout', action='store_false', help='disable VTKOutput')
    parser.set_defaults(vtkout=False)
    parser.add_argument('--plot', dest='plot', action='store_true', help='enable plotting')
    parser.add_argument('--no-plot', dest='plot', action='store_false', help='disable plotting')
    parser.set_defaults(plot=True)
    args = parser.parse_args()
    options = vars(args)
    print("call arguments: ", options)
    LaplaceBeltramiOnSphere(reflvls=options['reflvls'], vtkout=options['vtkout'], plot=options['plot'])    
        
