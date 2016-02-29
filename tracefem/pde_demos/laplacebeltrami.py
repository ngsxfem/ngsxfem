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


from SPDEutils import *
from xfem.lsetcurv import *


def InitialSolve(data):
    with TaskManager():
        SolveProblem(data)
    PostProcess(data)
    
def RefineAndSolve(data):
    problemdata = data["problemdata"]
    results = data["results"]
    
    with TaskManager():
        problemdata["Mesh"].Refine()
    results["level"] = results["level"] + 1
    with TaskManager():
        SolveProblem(data)
    PostProcess(data)

def LaplaceBeltramiOnSphere(reflvls = 1):
    options = MakeDefaultOptions();
    problemdata = Make3DProblem();
    discparams = MakeDefaultDiscretizationParams();
    pdeobj = MakePDEObjects(problemdata["Mesh"], problemdata, discparams);
    lsetmeshadap = LevelSetMeshAdaptation(problemdata["Mesh"], order=discparams["Order"], threshold=1000,
                                      discontinuous_qn=discparams["DiscontinuousQuasiNormal"])
    if not os.path.exists("LapBeltrResults"):
        os.makedirs("LapBeltrResults")

    results = MakeResultContainer(pdeobj, lsetmeshadap, problemdata, options);

    data = {"options" : options, "problemdata" : problemdata, "discparams" : discparams, "pdeobj" : pdeobj, "lsetmeshadap" : lsetmeshadap, "results" : results}


    
    InitialSolve(data)
    for i in range(reflvls-1):
        RefineAndSolve(data)
    if (reflvls > 1):
        PlotConvergence(data)
    
if __name__ == "__main__":
    if (len(sys.argv) > 1):
        LaplaceBeltramiOnSphere(reflvls = int(sys.argv[1]))
    else:
        LaplaceBeltramiOnSphere(reflvls = 3)
        
