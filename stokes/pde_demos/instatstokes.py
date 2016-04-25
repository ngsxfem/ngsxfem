from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.ngstd import *
from ngsolve.la import *

#import numpy as np

from math import sqrt

from time import sleep

from xfem.basics import *
from libngsxfem_xfem import *
from libngsxfem_xstokes import *
#from libngsxfem_levelset import *

material_params_artificial ={"eta" : [1,1], 
                             "rho" : [0.1,0], 
                             "g" : 9.81, 
                             "sigma" : 0.1}

material_params_water_butanol = {"eta" : [0.001388,0.003281], 
                                 "rho" : [986.506,845.442], 
                                 "g" : 9.81, 
                                 "sigma" : 0.00163}

solver_params_water_butanol_default ={ "lsetorder" : 1,
                                       "velorder" : 2, 
                                       "ref_space" : 1, 
                                       "dt" : 0.0001, 
                                       "timesteps" : 200,
                                       "relvelvec" : [0.0, -0.65, 0.0],
                                       "gradrepair" : False,
                                       "initial_cond" : "stokes"
                                     }

bubble_shape_2D_radius_0_06 = "( sqrt((x-0.0)*(x-0.0)+(y-0.0)*(y-0.0)) - 0.0036)"
bubble_shape_3D_radius_0_06 = "( sqrt((x-0.0)*(x-0.0)+(y-0.0)*(y-0.0)+(z-0.0)*(z-0.0)) - 0.0036)"

# from libxfem import *
class Levelset:
    def __init__(self, mesh, init_lset, order=1):
        self.fes = FESpace ("HDG", mesh, order=order,
                            flags = {"dirichlet" : [1,2,3,4]} )
        self.fes.Update()
        self.gf = GridFunction(self.fes,"levelset(disc.)")
        self.gf.Update()
        self.gf.components[0].Set(init_lset)
        self.gf.components[1].Set(init_lset, boundary=True)
        self.coef = GFCoeff(self.gf)

        self.cont_fes = FESpace ("h1ho", mesh, order=order)
        self.cont_fes.Update()
        self.cont_gf = GridFunction(self.cont_fes,"levelset(cont.)")
        self.cont_gf.Update()
        self.cont_gf.Set(init_lset)
        self.cont_coef = GFCoeff(self.cont_gf)

class VelocityPressure:
    def __init__(self, mesh, lset,
                 init_vals=[ConstantCF(0), ConstantCF(0)], 
                 dirichlet_vel = [1,2,3,4],
                 ref_space = 1,
                 order=1):
        self.fes = FESpace ("xstokes", mesh=mesh, order=order, 
                            flags = {"dirichlet_vel" : dirichlet_vel, "empty_vel" : True, 
                                     "dgjumps" : False, "ref_space" : ref_space})

        CastToXStokesFESpace(self.fes).SetLevelSet(lset.cont_coef)
        self.fes.Update()

        self.gf = GridFunction(self.fes,"velocity-pressure")
        self.gf.Update()
        # self.gf.components[0].Visualize("velocity(x)")
        # self.gf.components[1].Visualize("velocity(y)")
        # self.gf.components[2].Visualize("pressure")

        self.gf.components[0].components[0].Set(coefficient = init_vals[0], boundary = False)
        self.gf.components[1].components[0].Set(coefficient = init_vals[1], boundary = False)

        # self.gf.components[0].components[0].Set(coefficient = dirichlet_vel_vals[0], boundary = True)
        # self.gf.components[1].components[0].Set(coefficient = dirichlet_vel_vals[1], boundary = True)

        self.coef = [GFCoeff(self.gf.components[i]) for i in range(0,spacedim)]
    
class LevelsetSolver:
    def __init__(self, lset, velpre, dt, relvelvec = [0.0,0.0,0.0], gradrepair = False):
        self.dt = dt
        self.lset = lset 
        self.gradrepair = gradrepair
        self.f = LinearForm (lset.fes)
        lfi_inner = LFI (name = "source", dim = spacedim,  coef = [GFCoeff(lset.gf)])
        lfi_block = CompoundLFI ( lfi = lfi_inner, comp = 0 )
        self.f.Add (lfi_block)

        if (self.gradrepair):
            self.g = LinearForm (lset.fes)
            if (spacedim == 2):
                self.g.Add (LFI (name = "lsetcorr", dim = spacedim, coef = [GFCoeff(lset.gf),velpre.coef[0],velpre.coef[1]]))
            else:
                self.g.Add (LFI (name = "lsetcorr", dim = spacedim, coef = [GFCoeff(lset.gf),velpre.coef[0],velpre.coef[1],velpre.coef[2]]))

        self.a = BilinearForm (lset.fes, flags = { }, symmetric = False)
        bfi_inner = BFI (name = "mass", dim = spacedim,  coef = [ConstantCF(1.0/self.dt)])
        bfi_block = CompoundBFI (bfi = bfi_inner, comp = 0 )
        self.a.Add (bfi_block)

        self.a.Add (BFI (name = "HDG_convection", dim = spacedim, coef = velpre.coef))
        if (relvelvec[0] != 0 or relvelvec[1] != 0 or relvelvec[2] != 0):
            self.relvel = True
            self.a.Add (BFI (name = "HDG_convection", dim = spacedim, 
                             coef = [ConstantCF(relvelvec[i]) for i in range(0,spacedim)]))
        else:
            self.relvel = False

      
        self.c = Preconditioner (self.a, "direct", { "inverse" : "pardiso" })

        # lset project: definition of (bi-) and linear forms

        self.proj_f = LinearForm (lset.cont_fes)
        self.proj_f.Add (LFI (name = "source", dim = spacedim,  coef = [GFCoeff(lset.gf)]))

        self.proj_a = BilinearForm (lset.cont_fes, flags = { "symmetric" : True })
        self.proj_a.Add (BFI (name = "mass", dim = spacedim,  coef = [ConstantCF(1)]))

        self.proj_c = Preconditioner (self.proj_a, "direct", { "inverse" : "pardiso" })

        self.proj_a.Assemble(heapsize=10000000)
        self.proj_c.Update()

    def Update(self):
        self.f.Assemble(heapsize=10000000)
        if (self.gradrepair):
            self.g.Assemble(heapsize=10000000)
        self.a.Assemble(heapsize=10000000)

    def Solve(self, printrates = False, printrates_proj = False):
        self.c.Update()
        solver = CGSolver (self.a.mat, self.c.mat, printrates=printrates)

        if (self.relvel):
            tmp = self.lset.gf.vec.CreateVector()
            tmp2 = self.lset.gf.vec.CreateVector()
            tmp.data = 1.0/self.dt * self.f.vec.data
            if (self.gradrepair):
                tmp.data += self.g.vec.data
            tmp.data = tmp.data - self.a.mat *  self.lset.gf.vec.data
            tmp2.data = 0.0 * tmp2.data
            tmp2.data = solver * tmp.data
            self.lset.gf.vec.data = self.lset.gf.vec.data + tmp2.data
        else:
            if (self.gradrepair):
                rhs = self.lset.gf.vec.CreateVector()
                rhs.data = ((1.0/self.dt * self.f.vec.data) + self.g.vec.data)
                self.lset.gf.vec.data = solver * rhs.data
            else:
                self.lset.gf.vec.data = 1.0/self.dt * solver * self.f.vec

        self.proj_f.Assemble(heapsize=10000000)
        solver_proj = CGSolver (self.proj_a.mat, self.proj_c.mat, printrates=printrates_proj)
        self.lset.cont_gf.vec.data = solver_proj * self.proj_f.vec
        
class StokesSolver:
    def __init__(self, velpre, lset, dt,
                 params = {"eta" : [1,1], "rho" : [0.1,0], "g" : 9.81, "sigma" : 0.1},
                 initial = False
                 ):
        self.initial = initial
        self.dt = dt
        self.velpre = velpre
        self.f = LinearForm (self.velpre.fes)
        lfi_grav_inner = LFI (name = "xsource", dim = spacedim, coef = 
                              [ConstantCF(params["rho"][i]*params["g"]) for i in range(0,2)])
        self.f.Add (CompoundLFI ( lfi = lfi_grav_inner, comp=spacedim-1))
        self.f.Add (LFI (name = "xmodLBmeancurv", dim = spacedim, 
                         coef = [ConstantCF(params["sigma"]), lset.coef]))

        self.a = BilinearForm (self.velpre.fes, flags = { "symmetric" : True })
        self.a.Add (BFI (name = "xstokes", dim = spacedim, coef = [ConstantCF(params["eta"][0]), 
                                                                   ConstantCF(params["eta"][1])]))
        self.c = Preconditioner (self.a, "direct", { "inverse" : "pardiso" })


        if (self.initial == False):
            self.f_old = LinearForm (self.velpre.fes)
            for i in range(0,spacedim):
                lfi_inner_xi = LFI (name = "xsource", dim = spacedim,  
                                    coef = [self.velpre.coef[i],self.velpre.coef[i]])
                lfi_block_xi = CompoundLFI (lfi = lfi_inner_xi, comp = i )
                self.f_old.Add (lfi_block_xi)
                
                bfi_inner_xi = BFI (name = "xmass", dim = spacedim,  
                                    coef = [ConstantCF(1.0/self.dt),ConstantCF(1.0/self.dt)])
                bfi_block_xi = CompoundBFI (bfi = bfi_inner_xi, comp = i )
                self.a.Add (bfi_block_xi)

    def Update(self):
        self.velpre.fes.Update()
        # the order of the following operations is important as ....gf.Update() empties the vector if dimension has changed...
        if (self.initial == False): 
            self.f_old.Assemble(heapsize=10000000)
        self.velpre.gf.Update()
        self.f.Assemble(heapsize=10000000)
        self.a.Assemble(heapsize=10000000,reallocate=True)

    def Solve(self,printrates = False):
        self.c.Update()
        solver = CGSolver (self.a.mat, self.c.mat, printrates=printrates)

        if (self.initial == False):
            self.f.vec.data = self.f.vec.data + 1.0/self.dt * self.f_old.vec.data

        if (False): # non-homogenuous boundary conditions for the velocity...
            self.velpre.gf.components[1].components[0].Set(coefficient = ConstantCF(-0.1), boundary = True)
            self.f.vec.data = self.f.vec.data - self.a.mat *  self.velpre.gf.vec.data
            tmp = self.velpre.gf.vec.CreateVector()
            tmp = self.velpre.gf.vec
            self.velpre.gf.vec.data = solver * self.f.vec.data
            self.velpre.gf.vec.data = self.velpre.gf.vec.data + tmp
        else:
            self.velpre.gf.vec.data = solver * self.f.vec.data
        
       
def DoInstatStokes():
    print ("solve instationary stokes problem")
    
    matparams = material_params_water_butanol
    solverparams = solver_params_water_butanol_default

    global spacedim 
    spacedim = 2 # get that from mesh

    if (spacedim == 2):
        mesh = Mesh("d9_stokes.vol.gz")
        coef_initial_lset = VariableCF(bubble_shape_2D_radius_0_06)
    else:
        mesh = Mesh("d9_stokes_3D.vol.gz")
        coef_initial_lset = VariableCF(bubble_shape_3D_radius_0_06)

    lset = Levelset(mesh = mesh, init_lset = coef_initial_lset, 
                    order = solverparams["lsetorder"])
    global velpre
    velpre = VelocityPressure(mesh, lset, dirichlet_vel=[1,2,3,4], 
                              ref_space=solverparams["ref_space"], 
                              order=solverparams["velorder"]-1)

    lset_solver = LevelsetSolver(lset,velpre,dt=solverparams["dt"],
                                 relvelvec=solverparams["relvelvec"],
                                 gradrepair=solverparams["gradrepair"])

    if (solverparams["initial_cond"] == "stokes"):
        init_stokes_solver = StokesSolver(velpre,lset,params=matparams,dt=solverparams["dt"], initial = True)

        init_stokes_solver.Update()
        init_stokes_solver.Solve()

        print("after initial stokes solve")

    stokes_solver = StokesSolver(velpre,lset,params=matparams,dt=solverparams["dt"])

    print("before time loop - waiting for input(press enter) to start..")
    input()

    Draw (lset_solver.lset.gf, mesh, "lset" ) 
    Draw (CoefficientFunction(velpre.gf.components[0:2]), mesh, "vel" ) 
    
    for i in range(1,solverparams["timesteps"]+1):
        lset_solver.Update()
        lset_solver.Solve()
        stokes_solver.Update()
        stokes_solver.Solve()

        Redraw(blocking=True)
    
        print ('\rtime step ' + repr(i) + ' / ' + repr(solverparams["timesteps"]) + ' done. (ndof: ' + repr(velpre.fes.ndof) +')', end='');
    print("simulation finished")
    

DoInstatStokes()
