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
    square.AddRectangle([-1,-1],[1,1],bc="dirbound")
    mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

    mu1 = 1.0
    mu2 = 10.0

    R = 2.0/3.0
    aneg = 1.0/mu1
    apos = 1.0/mu2 + (1.0/mu1 - 1.0/mu2)*exp(x*x+y*y-R*R)
    gammaf = 0.5
    q = gammaf - pi*R*R/4.0*gammaf
    problem = {"ViscosityInner" : mu1,
               "ViscosityOuter" : mu2,
               "GammaF" : gammaf,
               "Levelset" : sqrt(x*x+y*y)-R,
               "SourceInnerX" : (exp(-1*(x*x+y*y))*((-8*y)+(4*x*x*y)+(4*y*y*y))+3*x*x),
               "SourceOuterX" : (exp(-1*(x*x+y*y))*((-8*y)+(4*x*x*y)+(4*y*y*y))+3*x*x),
               "SourceInnerY" : (exp(-1*(x*x+y*y))*((-4*x*x*x)+(8*x)-(4*x*y*y))),
               "SourceOuterY" : (exp(-1*(x*x+y*y))*((-4*x*x*x)+(8*x)-(4*x*y*y))),
               "SolutionOuterVelX" : (apos*exp(-1.0 *( x * x + y * y)) * -1.0 * y),
               "SolutionOuterVelY" : (apos*exp(-1.0 *( x * x + y * y)) * x),
               "SolutionInnerVelX" : (aneg*exp(-1.0 *( x * x + y * y)) * -1.0 * y),
               "SolutionInnerVelY" : (aneg*exp(-1.0 *( x * x + y * y)) * x),
               "SolutionInnerPressure" : (x*x*x + q),
               "SolutionOuterPressure" : (x*x*x - (pi*R*R/4.0*gammaf)),
               "NitscheParam" : 20,
               "GhostPenaltyParam" : -0.1,
               "Mesh" : mesh
              }
    return problem;

problemdata = Make2DProblem()
mesh = problemdata["Mesh"]
order = 1 # Pk+1 Pk

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order+1, threshold=1000, discontinuous_qn=True)

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

# c = Preconditioner(a, type="direct", flags= { "inverse" : "pardiso" })

uvp = GridFunction(Vh)

def SolveProblem():
    # Calculation of the deformation:
    deformation = lsetmeshadap.CalcDeformation(problemdata["Levelset"])
    # Applying the mesh deformation
    mesh.SetDeformation(deformation)

    Vh.Update()
    print("XStokesFESpace NDof:", Vh.ndof)
    uvp.Update()

    uvp.components[0].components[0].Set(problemdata["SolutionOuterVelX"], boundary=True, definedon=mesh.Boundaries("dirbound"))
    uvp.components[1].components[0].Set(problemdata["SolutionOuterVelY"], boundary=True, definedon=mesh.Boundaries("dirbound"))
    
    a.Assemble();
    f.Assemble();
    # c.Update();

    # the boundary value problem to be solved on each level
    inv = a.mat.Inverse(Vh.FreeDofs(), inverse="pardiso")
    f.vec.data -= a.mat * uvp.vec
    uvp.vec.data += inv * f.vec;
   
    sol_velocity_x = IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterVelX"],problemdata["SolutionInnerVelX"])
    sol_velocity_y = IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterVelY"],problemdata["SolutionInnerVelY"])
    sol_velocity = CoefficientFunction((sol_velocity_x,sol_velocity_y))
    
    velocity = CoefficientFunction (uvp.components[0:2])

    err_velocity_vec = velocity - sol_velocity
    err_velocity_sqr = err_velocity_vec * err_velocity_vec

    l2error_vel = sqrt(IntegrateOnWholeDomain(lsetmeshadap.lset_p1,mesh,
                                              coef=err_velocity_sqr,
                                              order=2*order+2))
    
    print("(velocity) l2 error = {}".format(l2error_vel))



    pressure_offset = IntegrateOnWholeDomain(lsetmeshadap.lset_p1,mesh,coef=CoefficientFunction (uvp.components[2]),order=order) / 4.0
    print (" pressure offset is {}".format(pressure_offset))
    pressure = CoefficientFunction (uvp.components[2]) - pressure_offset

    sol_pressure = IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterPressure"],problemdata["SolutionInnerPressure"])
    
    err_pressure_sqr = (pressure - sol_pressure)*(pressure - sol_pressure)

    l2error_pre = sqrt(IntegrateOnWholeDomain(problemdata["Levelset"],mesh,
                                              coef=err_pressure_sqr,
                                              order=2*order+1))
    print("(pressure) l2 error = {}".format(l2error_pre))


    Draw(problemdata["Levelset"],mesh,"levelset")
    Draw(lsetmeshadap.deform,mesh,"deformation")
    Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")
    Draw(velocity,mesh,"velocity")
    Draw(pressure,mesh,"pressure")

    mesh.UnsetDeformation()

def RefineAndSolve(n=1):
    for i in range(n):
        mesh.Refine()
        SolveProblem()
    
if __name__ == "__main__":
    SolveProblem()
    RefineAndSolve(2)
