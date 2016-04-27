from math import pi
# ngsolve stuff
from ngsolve import *

from xfem.utils import *

h = specialcf.mesh_size

# 3D: circle configuration
def Make3DProblem_Diffusion():
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    cube = CSGeometry()
    cube.Add (OrthoBrick(Pnt(-1.41,-1.41,-1.41), Pnt(1.41,1.41,1.41)))
    # mesh = Mesh (cube.GenerateMesh(maxh=0.5, quad_dominated=False))
    mesh = Mesh (cube.GenerateMesh(maxh=1, quad_dominated=False))
    mesh.Refine()
    
    a = 1
    c = 1
    
    problem = {"Diffusion" : a,
               "Convection" : None,
               "Reaction" : c,
               "Source" : (sin(pi*z)*(a*pi*pi*(1-z*z)+c)+a*cos(pi*z)*2*pi*z),
               "Solution" : sin(pi*z),
               "GradSolution" : CoefficientFunction((pi*cos(pi*z)*(-x*z),pi*cos(pi*z)*(-y*z),pi*cos(pi*z)*(1-z*z))),
               "VolumeStabilization" : a/h+c*h,
               "Levelset" : LevelsetExamples["sphere"],
               "GradLevelset" : CoefficientFunction((x,y,z)),
               "Lambda" : 10,
               "Iterative" : False,
               "Order" : 2,
               "Mesh" : mesh,
               "StaticCondensation" : True,
               "HDG": True,
               # "checkDGpattern" : True,
               # "checkCGpattern" : True,
               # "checkCGGPpattern" : True,
               # "checkHDGpattern" : True,
    }
    return problem;

def Make3DProblem_Convection():
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    cube = CSGeometry()
    cube.Add (OrthoBrick(Pnt(-1.41,-1.41,-1.41), Pnt(1.41,1.41,1.41)))
    # mesh = Mesh (cube.GenerateMesh(maxh=0.5, quad_dominated=False))
    mesh = Mesh (cube.GenerateMesh(maxh=1, quad_dominated=False))
    mesh.Refine()

    eps = 0.025
    c = 1
    problem = {"Diffusion" : 0.0,
               "Convection" : CoefficientFunction((-y*sqrt(1-z*z),x*sqrt(1-z*z),0)),
               "Reaction" : c,
               # "Source" : (12 * eps*sqrt(eps)*x*y*z/(eps+4*y*y) + 16 *eps*sqrt(eps)*(1-z*z)*x*y*z/((eps+4*z*z)*(eps+4*z*z)) + (6*eps*x*y+sqrt(x*x+y*y)*(x*x-y*y)+x*y)*atan(2*z/sqrt(eps))),
               "Source" : (sqrt(x*x+y*y)*(x*x-y*y)+x*y)*atan(2*z/sqrt(eps)) ,
               "Solution" : x*y*atan(2*z/sqrt(eps)),
               "GradSolution" : None,
               "VolumeStabilization" : 1+eps/h+c*h,
               "Levelset" : LevelsetExamples["sphere"],
               "GradLevelset" : CoefficientFunction((x,y,z)),
               "Lambda" : 10,
               "Iterative" : False,
               "Order" : 2,
               "Mesh" : mesh,
               "StaticCondensation" : True,
               "HDG": True,
               # "checkDGpattern" : True,
               # "checkHDGpattern" : True,
               # "checkCGpattern" : True,
               # "checkCGGPpattern" : True
    }
    return problem;


def Make3DProblem():
    return Make3DProblem_Diffusion()
    # return Make3DProblem_Convection()
