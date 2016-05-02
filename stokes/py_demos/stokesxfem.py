# solve the Laplace Beltrami equation on a sphere
# with Dirichlet boundary condition u = 0
from math import pi
# ngsolve stuff
from ngsolve import *
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

h = specialcf.mesh_size

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
               "SolutionOuterVelX_DX" : (2.0*x*y/mu2*exp(-1.0 *( x * x + y * y))),
               "SolutionOuterVelX_DY" : ((2.0*y*y/mu2-apos)*exp(-1.0 *( x * x + y * y))),
               "SolutionOuterVelY_DX" : ((apos - 2*x*x / mu2)*exp(-1.0 *( x * x + y * y))),
               "SolutionOuterVelY_DY" : (-2.0*x*y/mu2*exp(-1.0 *( x * x + y * y))),
               "SolutionInnerVelX_DX" : (aneg*2*x*y*exp(-1.0 *( x * x + y * y))),
               "SolutionInnerVelX_DY" : (aneg*(2*y*y-1)*exp(-1.0 *( x * x + y * y))),
               "SolutionInnerVelY_DX" : (aneg*(1-2*x*x)*exp(-1.0 *( x * x + y * y))),
               "SolutionInnerVelY_DY" : (aneg*(-2)*x*y*exp(-1.0 *( x * x + y * y))),
               "SolutionInnerPressure" : (x*x*x + q),
               "SolutionOuterPressure" : (x*x*x - (pi*R*R/4.0*gammaf)),
               "NitscheParam" : 20,
               "GhostPenaltyParam" : -0.1 * h * h, # GP-integrator has scaling h, but we need h^3
               "EmptyVel" : False,
               "Mesh" : mesh
              }
    return problem;

problemdata = Make2DProblem()
mesh = problemdata["Mesh"]
order = 1 # Pk+1 Pk

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order+1, threshold=1000, discontinuous_qn=True)

# ### Setting up discrete variational problem

dgjumps = "GhostPenaltyParam" in problemdata and problemdata["GhostPenaltyParam"] != 0.0
    
Vh = XStokesFESpace(mesh, order=order, levelset=lsetmeshadap.lset_p1, dirichlet=[1,2,3,4], dgjumps=dgjumps, empty_vel=problemdata["EmptyVel"])

a = BilinearForm(Vh, symmetric = True, flags = { })
a += TwoDomainStokesIntegrator(problemdata["ViscosityInner"],problemdata["ViscosityOuter"])
nitsche_a, nitsche_f = NitscheStokesIntegrators(problemdata["ViscosityInner"],
                                                problemdata["ViscosityOuter"],
                                                lamb=problemdata["NitscheParam"],
                                                gammaf=problemdata["GammaF"])

if dgjumps:
    a.components[2] += GhostPenaltyIntegrator(coefneg=1.0/problemdata["ViscosityInner"],
                                              coefpos=1.0/problemdata["ViscosityOuter"],
                                              stab_param=problemdata["GhostPenaltyParam"])
a += nitsche_a

f = LinearForm(Vh)
f += nitsche_f
f.components[0] += TwoDomainSourceIntegrator(problemdata["SourceInnerX"],
                                             problemdata["SourceOuterX"])
f.components[1] += TwoDomainSourceIntegrator(problemdata["SourceInnerY"],
                                             problemdata["SourceOuterY"])

# c = Preconditioner(a, type="direct", flags= { "inverse" : "pardiso" })

uvp = GridFunction(Vh)

def ApplyMeshTrafo():
    # Calculation of the deformation:
    deformation = lsetmeshadap.CalcDeformation(problemdata["Levelset"])
    # Applying the mesh deformation
    mesh.SetDeformation(deformation)

def UnapplyMeshTrafo():
    mesh.UnsetDeformation()
    
def SolveProblem():
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

def DictionaryAppend(cont,key,val):
    if (cont != None):
        if (not(key in cont)):
            cont[key] = []
        cont[key].append(val)

    
def ComputeErrors(statistics_dict=None):
    DictionaryAppend(statistics_dict,"total ndof",Vh.ndof)
    ### velocity error
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
    DictionaryAppend(statistics_dict,"error uv (L2)",l2error_vel)

    velocity_grad = CoefficientFunction ((uvp.components[0].Deriv()[0],uvp.components[0].Deriv()[1],uvp.components[1].Deriv()[0],uvp.components[1].Deriv()[1]))

    sol_velocity_grad = CoefficientFunction((
        IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterVelX_DX"],problemdata["SolutionInnerVelX_DX"]),
        IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterVelX_DY"],problemdata["SolutionInnerVelX_DY"]),
        IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterVelY_DX"],problemdata["SolutionInnerVelY_DX"]),
        IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterVelY_DY"],problemdata["SolutionInnerVelY_DY"])))

    err_velocity_grad_vec = sol_velocity_grad - velocity_grad
    err_velocity_grad_sqr = err_velocity_grad_vec * err_velocity_grad_vec

    h1seminormerror_vel = sqrt(IntegrateOnWholeDomain(lsetmeshadap.lset_p1,mesh,
                                              coef=err_velocity_grad_sqr,
                                              order=2*order))
    
    print("(velocity) h1 error = {}".format(l2error_vel + h1seminormerror_vel))
    DictionaryAppend(statistics_dict,"error uv (H1)",l2error_vel + h1seminormerror_vel)

        
    ### pressure offset correction
    pressure_offset = IntegrateOnWholeDomain(lsetmeshadap.lset_p1,mesh,coef=CoefficientFunction (uvp.components[2]),order=order) / 4.0
    print (" pressure offset is {}".format(pressure_offset))
    pressure = CoefficientFunction (uvp.components[2]) - pressure_offset

    ### pressure error
    sol_pressure = IfPos(lsetmeshadap.lset_p1,problemdata["SolutionOuterPressure"],problemdata["SolutionInnerPressure"])
    
    err_pressure_sqr = (pressure - sol_pressure)*(pressure - sol_pressure)

    l2error_pre = sqrt(IntegrateOnWholeDomain(problemdata["Levelset"],mesh,
                                              coef=err_pressure_sqr,
                                              order=2*order+1))
    print("(pressure) l2 error = {}".format(l2error_pre))
    DictionaryAppend(statistics_dict,"error p (L2)",l2error_pre)

def PrintDictionaryFormatted(statistics_dict):
    if statistics_dict == None:
        return
    else:
        for key,val in statistics_dict.items():
            nlvl = len(statistics_dict[key])
            break
    for key,val in statistics_dict.items():
        print(" {:>14} ".format(key),end="")
        if ("error" in key):
            print("( eoc) ".format(key),end="")
    print("")
    for key,val in statistics_dict.items():
        print("----------------",end="")
        if ("error" in key):
            print("-------".format(key),end="")
    print("")
    for lvl in range(0,nlvl):
        for key,val in statistics_dict.items():
            if ("error" in key):
                print("   {:10e} ".format(val[lvl]),end="")
                if (lvl>0):
                    print("({:.2f}) ".format(log(val[lvl-1]/val[lvl])/log(2)),end="")
                else:
                    print("( -- ) ",end="")
            else:
                print("    {:>11} ".format(val[lvl]),end="")
        print("")

    
def MakeVTKOutput(lvl=0):

    Vhh1 = H1(mesh, order=order+1, dirichlet=[1,2,3,4])
    Vhp = H1(mesh, order=order, dirichlet=[])

    Vhuvnegpos = FESpace([Vhh1,Vhh1])
    Vhpnegpos = FESpace([Vhp,Vhp])
    
    unegpos = GridFunction(Vhuvnegpos)
    vnegpos = GridFunction(Vhuvnegpos)
    pnegpos = GridFunction(Vhpnegpos)
    
    XToNegPos(uvp.components[0],unegpos)
    XToNegPos(uvp.components[1],vnegpos)
    XToNegPos(uvp.components[2],pnegpos)

    velneg = CoefficientFunction((unegpos.components[0],vnegpos.components[0]))
    velpos = CoefficientFunction((unegpos.components[1],vnegpos.components[1]))
    
    
    pairs = [(problemdata["Levelset"],"levelset"),
             (lsetmeshadap.deform,"deformation"),
             (lsetmeshadap.lset_p1,"levelset(P1)"),
             (velneg,"velocity_neg"),
             (velpos,"velocity_pos"),
             (pnegpos.components[0],"pressure_neg"),
             (pnegpos.components[1],"pressure_pos")]

    coefs,names = [list(t) for t in zip(*pairs)]
    VTKOutput(ma=mesh,coefs=coefs,names=names,filename="vtkout_lvl"+str(lvl),subdivision=3).Do()
    
        
def DrawSolution():
    Draw(problemdata["Levelset"],mesh,"levelset")
    Draw(lsetmeshadap.deform,mesh,"deformation")
    Draw(lsetmeshadap.lset_p1,mesh,"levelset(P1)")
    Draw(CoefficientFunction (uvp.components[0:2]),mesh,"velocity")
    Draw(CoefficientFunction (uvp.components[2]),mesh,"pressure")

def RefineAndSolve(n=1, statistics_dict=None):
    for i in range(n):
        if (i!=0):
            mesh.Refine()
        ApplyMeshTrafo()
        SolveProblem()
        ComputeErrors(statistics_dict)
        UnapplyMeshTrafo()
        MakeVTKOutput(lvl=i)

if __name__ == "__main__":
    statistics_dict=dict()
    RefineAndSolve(3,statistics_dict)
    PrintDictionaryFormatted(statistics_dict)
    DrawSolution()
