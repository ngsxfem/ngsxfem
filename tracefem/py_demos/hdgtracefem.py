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
# For Integration on Interface
from xfem.basics import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *
# For TraceFEM-Integrators (convenience)
from xfem.tracefem import *
# For HDGTraceFEM-Integrators (convenience)
from xfem.hdgtracefem import *

# 3D: circle configuration
def Make3DProblem():
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    cube = CSGeometry()
    cube.Add (OrthoBrick(Pnt(-1.41,-1.41,-1.41), Pnt(1.41,1.41,1.41)))
    # mesh = Mesh (cube.GenerateMesh(maxh=0.5, quad_dominated=False))
    mesh = Mesh (cube.GenerateMesh(maxh=1, quad_dominated=False))
    mesh.Refine()

    problem = {"Diffusion" : 1.0,
               "Convection" : None,
               "Reaction" : 1.0,
               "Source" : sin(pi*z)*(pi*pi*(1-z*z)+1)+cos(pi*z)*2*pi*z,
               "Solution" : sin(pi*z),
               "GradSolution" : CoefficientFunction((pi*cos(pi*z)*(-x*z),pi*cos(pi*z)*(-y*z),pi*cos(pi*z)*(1-z*z))),
               "SurfGradSolution" : CoefficientFunction((pi*cos(pi*z)*(-x*z),pi*cos(pi*z)*(-y*z),pi*cos(pi*z)*(1-z*z))),
               "VolumeStabilization" : 1,
               "Levelset" : sqrt(x*x+y*y+z*z)-1,
               "GradLevelset" : CoefficientFunction((x,y,z)),
               "Lambda" : 10,
               "Iterative" : True,
               "Order" : 2,
               "Mesh" : mesh,
               "StaticCondensation" : False,
               "HDG": False
              }
    return problem;

class Discretization(object):
    """
    Discretization: definition of FESpaces, bi/linear forms, etc..
    """
    
    def __init__(self, problemdata):
        """
        """
        self.problemdata = problemdata
        self.mesh = problemdata["Mesh"]
        self.order = problemdata["Order"]
        if (problemdata["VolumeStabilization"] and problemdata["StaticCondensation"]):
            self.static_condensation = True
        else:
            self.static_condensation = False
            
        self.lsetmeshadap = LevelSetMeshAdaptation(self.mesh, order=self.order, threshold=1000, discontinuous_qn=True,heapsize=10000000)

        ### Setting up discrete variational problem
        if self.problemdata["HDG"]:
            self.Vh_l2 = L2(self.mesh, order=self.order, dirichlet=[])
            Vh_l2_tr = TraceFESpace(self.mesh, self.Vh_l2, problemdata["Levelset"])
            self.Vh_facet = FacetFESpace(self.mesh, order=self.order, dirichlet=[], flags = {"highest_order_dc" : False})
            Vh_facet_tr = TraceFESpace(self.mesh, self.Vh_facet, problemdata["Levelset"])
            self.Vh_tr = FESpace([Vh_l2_tr,Vh_facet_tr])
            
            self.a = BilinearForm(self.Vh_tr, symmetric = True, flags = {"eliminate_internal" : self.static_condensation})
            if (problemdata["Reaction"] != None):
                self.a.components[0] += TraceMass(problemdata["Reaction"])
            if (problemdata["Diffusion"] != None):
                self.a += HDGTraceLaplaceBeltrami(problemdata["Diffusion"],
                                                  param_IP_edge = problemdata["Lambda"],
                                                  param_normaldiffusion = problemdata["VolumeStabilization"],
                                                  param_IP_facet = problemdata["Lambda"])
            
            self.f = LinearForm(self.Vh_tr)
            if (problemdata["Source"] != None):
                self.f.components[0] += TraceSource(problemdata["Source"])
            if (not problemdata["VolumeStabilization"]):
                raise Exception(" cannot turn off VolumeStabilization with HDG ")
            if (problemdata["Convection"] != None):
                raise Exception(" No HDG convection yet ")
        else:
            self.Vh = H1(self.mesh, order=self.order, dirichlet=[])
            self.Vh_tr = TraceFESpace(self.mesh, self.Vh, problemdata["Levelset"])

            self.a = BilinearForm(self.Vh_tr, symmetric = True, flags = {"eliminate_internal" : self.static_condensation})
            if (problemdata["Reaction"] != None):
                self.a += TraceMass(problemdata["Reaction"])
            if (problemdata["Diffusion"] != None):
                self.a += TraceLaplaceBeltrami(problemdata["Diffusion"])
            if (problemdata["VolumeStabilization"]!=None):
                self.a += NormalLaplaceStabilization(problemdata["Diffusion"]*problemdata["VolumeStabilization"],
                                                     # problemdata["GradLevelset"])
                                                     self.lsetmeshadap.lset_p1.Deriv())
            if (problemdata["Convection"] != None):
                self.a += TraceConvection(problemdata["Convection"])

            self.f = LinearForm(self.Vh_tr)
            if (problemdata["Source"] != None):
                self.f += TraceSource(problemdata["Source"])
        if (self.problemdata["Iterative"]):
            self.c = Preconditioner(self.a, type="local", flags= { "test" : True })
        else:
            self.c = Preconditioner(self.a, type="direct", flags= { "test" : True })
            
        self.u = GridFunction(self.Vh_tr)
  
    def UpdateSpace(self):
        if self.problemdata["HDG"]:
            self.Vh_facet.Update()
            self.Vh_l2.Update()
            self.Vh_tr.Update()
        else:
            self.Vh.Update()
            self.Vh_tr.Update()

        # if self.static_condensation:
        #     for i in range(self.Vh_tr.ndof):
        #         if (self.Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF):
        #             self.Vh_tr.FreeDofs()[i] = 0

    def VolumeSolution(self):
        if self.problemdata["HDG"]:
            return self.u.components[0]
        else:
            return self.u
                    
    def SolveProblem(self):
        results = {}
        # Calculation of the deformation:
        deformation = discretization.lsetmeshadap.CalcDeformation(self.problemdata["Levelset"])
        # Applying the mesh deformation
        self.mesh.SetDeformation(deformation)
        
        self.UpdateSpace()
        
        if (self.problemdata["HDG"]):
            Vh_l2 = L2(self.mesh, order=self.order, dirichlet=[])
            Vh_l2_tr = TraceFESpace(self.mesh, self.Vh_l2, problemdata["Levelset"], dgjumps=True)
            results["dg_ndofs"] = Vh_l2_tr.ndof
            b = BilinearForm(Vh_l2_tr)
            b.Assemble()
            results["dg_nze"] = b.mat.AsVector().size
            
        
        
        results["total_ndofs"] = self.Vh_tr.ndof
        global_ndofs = self.Vh_tr.ndof
        for i in range(self.Vh_tr.ndof):
        #     print (str(i) + " : " + str(Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF))
            if (self.Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF):
                # self.Vh_tr.FreeDofs()[i] = 0
                global_ndofs = global_ndofs - 1
        results["global_ndofs"] = global_ndofs
        
        self.u.Update()
        self.a.Assemble(heapsize=10000000);
        self.f.Assemble();
        
        self.c.Update();
        
        solvea = CGSolver( mat=self.a.mat, pre=self.c.mat, complex=False, printrates=False, precision=1e-8, maxsteps=200000)
        # u.vec.data = ainv * f.vec
        
        if self.static_condensation:
            self.f.vec.data += self.a.harmonic_extension_trans * self.f.vec
        self.u.vec.data = solvea * self.f.vec;
        if self.static_condensation:
            self.u.vec.data += self.a.inner_solve * self.f.vec
            self.u.vec.data += self.a.harmonic_extension * self.u.vec
            
        print("nze: " + str(self.a.mat.AsVector().size))
        results["nze"] = self.a.mat.AsVector().size
        
        last_num_its = solvea.GetSteps()
        print("number of iterations: " + str(last_num_its))
        results["numits"] = last_num_its
        
        coef_error_sqr = (self.VolumeSolution() - self.problemdata["Solution"])*(self.VolumeSolution() - self.problemdata["Solution"])
        l2diff = sqrt(IntegrateOnInterface(self.lsetmeshadap.lset_p1,self.mesh,coef_error_sqr,order=2*self.order+2))
        print("l2diff = {}".format(l2diff))
        results["l2err"] = l2diff
        
        
        nhelp = self.lsetmeshadap.lset_p1.Deriv()
        n = 1.0/sqrt(nhelp*nhelp) * nhelp
        un = self.VolumeSolution().Deriv()*n
        coef_gradnormal_error = un*un
        h1normdiff = sqrt(IntegrateOnInterface(self.problemdata["Levelset"],self.mesh,coef_gradnormal_error,order=2*self.order))
        print("h1normerr = {}".format(h1normdiff))
        results["h1normerr"] = h1normdiff
        
        tanggrad = self.VolumeSolution().Deriv() - un * n
        coef_gradtang_error = (tanggrad - self.problemdata["GradSolution"])*(tanggrad - self.problemdata["GradSolution"])
        h1tangdiff = sqrt(IntegrateOnInterface(self.problemdata["Levelset"],self.mesh,coef_gradtang_error,order=2*self.order))
        print("h1tangerr = {}".format(h1tangdiff))
        results["h1tangerr"] = h1tangdiff
        
        maxdistlset = self.lsetmeshadap.CalcMaxDistance(self.problemdata["Levelset"]);
        print("maxdist = {}".format(maxdistlset))
        results["maxdist"] = maxdistlset
        
        self.mesh.UnsetDeformation()
        return results

def PrintTimers(substring="HDG"):            
    ### print hdg-intergrator timers
    hdgtimers = [a for a in Timers() if substring in a["name"]]
    hdgtimers = sorted(hdgtimers, key=lambda k: k["time"], reverse=True)
    for timer in hdgtimers:
        print("{:<45}: {:6} cts, {:.8f} s, {:.6e} s(avg.)"
              .format(timer["name"],timer["counts"],timer["time"],timer["time"]/timer["counts"]))

import argparse        
import pickle
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Solve Laplace-Betrami example problem.')
    parser.add_argument('--reflvls', metavar='N', type=int, default=2, help='number of refinement levels (>=1)')
    parser.add_argument('--vtkout', dest='vtkout', action='store_true', help='enable VTKOutput')
    parser.add_argument('--no-vtkout', dest='vtkout', action='store_false', help='disable VTKOutput')
    parser.set_defaults(vtkout=False)
    args = parser.parse_args()
    options = vars(args)
    print("call arguments: ", options)
    
    # if True:
    with TaskManager():
        problemdata = Make3DProblem()
        discretization = Discretization(problemdata)

        orders = [4]
        l2diffresults = []
        resultdict = {}
        for i in range(options['reflvls']):
            l2diffresults.append([])
            if (i!=0):
                discretization.mesh.Refine()
            for order in orders:
                print("""
                ------------------------------------------------------------------
                   refinement level = {}, order = {}
                ------------------------------------------------------------------
                """.format(i,order))
                
                problemdata["Order"] = order
                result = False
                while result == False:
                    try:
                        discretization = Discretization(problemdata)
                        resultdict[(i,order)] = discretization.SolveProblem()
                        print("simulation results:\n {}".format(resultdict[(i,order)]))
                        l2diffresults[i].append(resultdict[(i,order)]["l2err"])
                        result = True
                    except Exception as e:
                        print ("Exception = ".format(e))
    #                    l2diffresults[i].append(None)
                        result = False
                
                with open('results.pkl', 'wb') as f:
                    pickle.dump(resultdict, f, 0) #pickle.HIGHEST_PROTOCOL)
                print(l2diffresults)
                print(list(map(list, zip(*l2diffresults))))
            RefineAtLevelSet(gf=discretization.lsetmeshadap.lset_p1)
    print(list(map(list, zip(*l2diffresults))))


    for i in range(options['reflvls']):
        for order in orders:
            print("({},{}): {}".format(i,order,resultdict[(i,order)]))
            
    # print(resultdict)
    #Draw(discretization.u.components[0],discretization.mesh,"u",draw_surf=False)
    

    # with open('results.pkl', 'rb') as f:
    #     resultdict = pickle.load(f)

    PrintTimers(substring="HDG")
    PrintTimers(substring="TraceFEM")
    PrintTimers(substring="XFE")
    
    if (options['vtkout']):
        vtk = VTKOutput(ma=discretization.mesh,coefs=[discretization.lsetmeshadap.lset_p1,
                                       discretization.lsetmeshadap.deform,u.components[0]],
                        names=["lsetp1","deform","u"],filename="vtkout_",subdivision=1)
        vtk.Do()

    
#observations:
# 1. for p=1 it seems that lsetcurving p+1=2 is necessary to achieve optimal order convergence (not true for p>1)
# 2. l2order=p+1, facetorder=p gives stable result (optimal order, i.e. p+1 in the l2-norm) - however p+1 together with highest_order_dc is not stable yet (?!)


#L2diff
#p1-> [   0.8096913452313316,      0.30840257878176175,    0.11076925787923653,   0.059828714229247425],
#p2-> [  0.11632850331153247,     0.011519875896377813,   0.001350445617732032,  0.0001810449912049936],
#p3-> [  0.03555471657059834,    0.0028774186159388226, 0.00017200706248403185, 1.1569149801119469e-05],
#p4-> [  0.00998605191678209,   0.00018175774468649536,  7.737886033220925e-06, 3.3478017691898923e-07],
#p5-> [0.0025838731998405035,   2.2492819637261814e-05,  9.017758296511354e-07,                       ],

