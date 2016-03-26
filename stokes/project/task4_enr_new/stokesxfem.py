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
# basic xfem functionality
from xfem.basics import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *
# For Stokes-FESpace and Stokes-Integrators (convenience)
from xfem.stokes import *

# 2D: circle configuration
def Make2DProblem():
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([-1,-1],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=2, quad_dominated=False))
    mesh.Refine()
    mesh.Refine()

    problem = {"ViscosityInner" : 1.0,
               "ViscosityOuter" : 10.0,
               "GammaF" : 0.5,
               "Levelset" : sqrt(x*x+y*y)-2.0/3.0,
               "SourceInnerX" : (exp(-1*(x*x+y*y))*((-8*y)+(4*x*x*y)+(4*y*y*y))+3*x*x),
               "SourceOuterX" : (exp(-1*(x*x+y*y))*((-8*y)+(4*x*x*y)+(4*y*y*y))+3*x*x),
               "SourceInnerY" : (exp(-1*(x*x+y*y))*((-4*x*x*x)+(8*x)-(4*x*y*y))),
               "SourceOuterY" : (exp(-1*(x*x+y*y))*((-4*x*x*x)+(8*x)-(4*x*y*y))),
               "NitscheParam" : 20,
               "GhostPenaltyParam" : -0.1,
               "Mesh" : mesh
              }
    return problem;

problemdata = Make2DProblem()
mesh = problemdata["Mesh"]
order = 1 # Pk+1 Pk

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)

# ### Setting up discrete variational problem

Vh = XStokesFESpace(mesh, order=order, levelset=problemdata["Levelset"], dirichlet=[1,2,3,4])

a = BilinearForm(Vh, symmetric = True, flags = { })
a += TwoDomainStokesIntegrator(problemdata["ViscosityInner"],problemdata["ViscosityOuter"])
nitsche_a, nitsche_f = NitscheStokesIntegrators(problemdata["ViscosityInner"],
                                                problemdata["ViscosityOuter"],
                                                lamb=problemdata["NitscheParam"],
                                                gammaf=problemdata["GammaF"])
a += nitsche_a

f = LinearForm(Vh)
f += nitsche_f
f.components[0] += TwoDomainSourceIntegrator(problemdata["SourceInnerX"],
                                             problemdata["SourceOuterX"])
f.components[1] += TwoDomainSourceIntegrator(problemdata["SourceInnerY"],
                                             problemdata["SourceOuterY"])

c = Preconditioner(a, type="direct", flags= { "inverse" : "pardiso" })

uvp = GridFunction(Vh)

def SolveProblem():
    # Calculation of the deformation:
    deformation = lsetmeshadap.CalcDeformation(problemdata["Levelset"])
    # Applying the mesh deformation
    mesh.SetDeformation(deformation)

    Vh.Update()
    print("XStokesFESpace NDof:", Vh.ndof)
    uvp.Update()

    u.components[0].Set(uin, definedon=mesh.Boundaries("inlet"))
    
    a.Assemble();
    f.Assemble();
    c.Update();

    inv = a.mat.Inverse(Vh.FreeDofs(), inverse="pardiso")
    # # the boundary value problem to be solved on each level
    uvp.vec.data = inv * f.vec;
   
#     mesh.UnsetDeformation()

if __name__ == "__main__":
    SolveProblem()
    velocity = CoefficientFunction (uvp.components[0:2])
    absvelocity = sqrt(velocity*velocity)
    Draw(velocity,mesh,"velocity")
    Draw(absvelocity,mesh,"absvelocity")
    print("Boundary conditions not yet set...")
    print("Error not yet evaluated...")
    
