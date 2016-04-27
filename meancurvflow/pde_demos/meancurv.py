from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.ngstd import *
from ngsolve.la import *

import numpy as np

# from math import sin
from time import sleep

from xfem.basics import *

# from libxfem import *

class npMeanCurvFlow(PyNumProc):

    def Do(self, heap):
        
        print ("solve mean curvature flow")
        
        print ("needs to be updated anyway....")

        input()
        #set time step
        dt = 0.5
        timesteps = 150

        lsetorder = 4
        pde = self.pde

        mesh = pde.Mesh()

        coef_initial_lset = pde.coefficients["initial_lset"]

        print(" initializing all fespaces, gridfunctions, etc..")

        # lset functions (fes,gf)

        # if (True): # initial lset problem
        fes_lset = FESpace ("HDG", mesh, order=lsetorder)
        fes_lset.Update()
        gf_lset = GridFunction(fes_lset)
        gf_lset.Update()
        gf_lset.component(1).Set(coef_initial_lset)

        fes_cont_lset = FESpace ("h1ho", mesh, order=lsetorder)
        fes_cont_lset.Update()
        gf_cont_lset = GridFunction(fes_cont_lset)
        gf_cont_lset.Update()
        gf_cont_lset.Set(coef_initial_lset)
        coef_lset = GFCoeff(gf_cont_lset)

        # stokes functions  (fes,gf)

        fescomp = FESpace ("xstokes", mesh, order=1, flags = {"dirichlet_vel" : [1,2,3,4], "empty_vel" : True, "dgjumps" : False, "ref_space" : 1})

        CastToXStokesFESpace(fescomp).SetLevelSet(coef_lset)
        fescomp.Update()

        uvp = GridFunction(fescomp)
        uvp.Update()

        # dirichlet boundary conditions for velocity ... 
        # uvp.component(2).component(1).Set(coefficient = ConstantCF(-0.1), boundary = True)

        coef_velx = GFCoeff(uvp.component(1))
        coef_vely = GFCoeff(uvp.component(2))


        # coef_velx = ConstantCF(1)
        # coef_vely = GFCoeff(uvp.component(2))


        print(" initializing all fespaces, gridfunctions, etc.. done \n")

        print(" defining all (bi)linear forms..")
        # lset transport: definition of (bi-) and linear forms


        f_lset = LinearForm (fes_lset)
        lfi_inner = LFI (name = "source", dim = 2,  coef = [GFCoeff(gf_lset)])
        lfi_block = LFI (type = "compound", lfi = lfi_inner, comp = 1 )
        f_lset.Add (lfi_block)
        
        a_lset = BilinearForm (fes_lset, flags = { "symmetric" : False })
        bfi_inner = BFI (name = "mass", dim = 2,  coef = [ConstantCF(1.0/dt)])
        bfi_block = BFI (type = "compound", bfi = bfi_inner, comp = 1 )
        a_lset.Add (bfi_block)

        a_lset.Add (BFI (name = "HDG_convection", dim = 2, coef = [coef_velx, coef_vely]))
        
        c_lset = Preconditioner (a_lset, "direct", { "inverse" : "pardiso" })
        


        coef_disc_lset = GFCoeff(gf_lset)

        # lset project: definition of (bi-) and linear forms

        f_cont_lset = LinearForm (fes_cont_lset)
        f_cont_lset.Add (LFI (name = "source", dim = 2,  coef = [coef_disc_lset]))
        f_cont_lset.Assemble()
        
        a_cont_lset = BilinearForm (fes_cont_lset, flags = { "symmetric" : True })
        a_cont_lset.Add (BFI (name = "mass", dim = 2,  coef = [ConstantCF(1)]))

        c_cont_lset = Preconditioner (a_cont_lset, "direct", { "inverse" : "pardiso" })


        # Stokes Kram

        f = LinearForm (fescomp)
        lfi_grav_inner = LFI (name = "xsource", dim = 2, coef = [ConstantCF(1), ConstantCF(0)])
        f.Add (LFI ( type = "compound", lfi = lfi_grav_inner, comp=2))
        f.Add (LFI (name = "xmodLBmeancurv", dim = 2, coef = [ConstantCF(0.1), coef_lset]))

        a = BilinearForm (fescomp, flags = { "symmetric" : True })
        a.Add (BFI (name = "xstokes", dim = 2, coef = [ConstantCF(1), ConstantCF(1)]))
        
        print(" defining all (bi)linear forms.. done\n")

        # input()

        #solve initial lset

        f_lset.Assemble()
        a_lset.Assemble()
        c_lset.Update()

        solver_lset = CGSolver (a_lset.mat, c_lset.mat, printrates=True)
        gf_lset.vec.data = 1.0/dt * solver_lset * f_lset.vec

        #project initial lset (-> continuous)

        a_cont_lset.Assemble()
        c_cont_lset.Update()

        solver_cont_lset = CGSolver (a_cont_lset.mat, c_cont_lset.mat, printrates=False)
        gf_cont_lset.vec.data = solver_cont_lset * f_cont_lset.vec

        #solve initial stokes

        f.Assemble()
        c = Preconditioner (a, "direct", { "inverse" : "pardiso" })
        a.Assemble()

        c.Update()
        solver = CGSolver (a.mat, c.mat, printrates=True)
        uvp.vec.data = solver * f.vec


        # a_cont_lset_vis = BilinearForm (fes_cont_lset, flags = {"nonassemble" : True })
        # bfi_inner_vis = BFI (name = "mass", dim = 2,  coef = [ConstantCF(1.0)])
        # a_cont_lset_vis.Add(bfi_inner_vis)

        # df = DrawFlux(bf = a_cont_lset_vis, gf = gf_cont_lset, label = "test-vis", applyd = True)
        # input()


        print("after initial stokes solve")

        f_old = LinearForm (fescomp)

        #todo add mass and source term to stokes
        lfi_inner_x = LFI (name = "xsource", dim = 2,  coef = [coef_velx,coef_velx])
        lfi_block_x = LFI (type = "compound", lfi = lfi_inner_x, comp = 1 )
        f_old.Add (lfi_block_x)

        lfi_inner_y = LFI (name = "xsource", dim = 2,  coef = [coef_vely,coef_vely])
        lfi_block_y = LFI (type = "compound", lfi = lfi_inner_y, comp = 2 )
        f_old.Add (lfi_block_y)

        bfi_inner_x = BFI (name = "xmass", dim = 2,  coef = [ConstantCF(1.0/dt),ConstantCF(1.0/dt)])
        bfi_block_x = BFI (type = "compound", bfi = bfi_inner_x, comp = 1 )
        a.Add (bfi_block_x)

        bfi_inner_y = BFI (name = "xmass", dim = 2,  coef = [ConstantCF(1.0/dt),ConstantCF(1.0/dt)])
        bfi_block_y = BFI (type = "compound", bfi = bfi_inner_y, comp = 2 )
        a.Add (bfi_block_y)


        input()
        for i in range(1,timesteps+1):
            
            f_lset.Assemble()
            a_lset.Assemble()
            c_lset.Update()
            
            # print("solve lset")
            solver_lset = CGSolver (a_lset.mat, c_lset.mat, printrates=False)
            gf_lset.vec.data = 1.0/dt * solver_lset * f_lset.vec
            # print("solve lset done\n")
            
            # print("project lset")
            f_cont_lset.Assemble()
            gf_cont_lset.vec.data = solver_cont_lset * f_cont_lset.vec
            # print("project lset done\n")

            Redraw(blocking=True)
        
            # sleep(1)

            fescomp.Update()
            uvp.Update()
            
            f.Assemble()
            f_old.Assemble()
            c = Preconditioner (a, "direct", { "inverse" : "pardiso" })
            a.Assemble(reallocate=True)
            
            c.Update()
            solver = CGSolver (a.mat, c.mat, printrates=False)
            f.vec.data = f.vec.data + 1.0/dt * f_old.vec.data
            uvp.vec.data = solver * f.vec

            Redraw(blocking=True)
        
            print ('\rtime step ' + repr(i) + ' / ' + repr(timesteps) + ' done.', end='');
            # sleep(1)

        input()


        
