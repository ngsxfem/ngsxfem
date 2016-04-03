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
    mesh = Mesh (cube.GenerateMesh(maxh=0.5, quad_dominated=False))
    # mesh.Refine()
    # mesh.Refine()
    # mesh.Refine()
    # mesh.Refine()
    # mesh.Refine()

    problem = {"Diffusion" : 1.0,
               "Convection" : None,
               "Reaction" : 1.0,
               "Source" : sin(pi*z)*(pi*pi*(1-z*z)+1)+cos(pi*z)*2*pi*z,
               "Solution" : sin(pi*z),
               "GradSolution" : CoefficientFunction((pi*cos(pi*z)*(-x*z),pi*cos(pi*z)*(-y*z),pi*cos(pi*z)*(1-z*z))),
               "SurfGradSolution" : CoefficientFunction((pi*cos(pi*z)*(-x*z),pi*cos(pi*z)*(-y*z),pi*cos(pi*z)*(1-z*z))),
               "VolumeStabilization" : 1,
               "Levelset" : sqrt(x*x+y*y+z*z)-1,
               "Lambda" : 10,
               "Order" : 2,
               "Mesh" : mesh
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
        if (problemdata["VolumeStabilization"]):
            self.static_condensation = True
        else:
            self.static_condensation = False
            
        self.lsetmeshadap = LevelSetMeshAdaptation(self.mesh, order=self.order, threshold=1000, discontinuous_qn=True,heapsize=10000000)

        ### Setting up discrete variational problem
        with TaskManager():
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
  
        self.c = Preconditioner(self.a, type="local", flags= { "test" : True })
  
        self.u = GridFunction(self.Vh_tr)
  

    def SolveProblem(self):
        with TaskManager():
        
            # Calculation of the deformation:
            deformation = discretization.lsetmeshadap.CalcDeformation(self.problemdata["Levelset"])
            # Applying the mesh deformation
            self.mesh.SetDeformation(deformation)
            
            self.Vh_facet.Update()
            self.Vh_l2.Update()
            self.Vh_tr.Update()
            
            for i in range(self.Vh_tr.ndof):
            #     print (str(i) + " : " + str(Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF))
                if (self.Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF):
                    self.Vh_tr.FreeDofs()[i] = 0
            
            self.u.Update()
            self.a.Assemble(heapsize=10000000);
            self.f.Assemble();
            
            self.c.Update();
            
            solvea = CGSolver( mat=self.a.mat, pre=self.c.mat, complex=False, printrates=False, precision=1e-9, maxsteps=200000)
            
            # solvea = self.a.mat.Inverse(self.Vh_tr.FreeDofs())
            # u.vec.data = ainv * f.vec
            
            if self.static_condensation:
                self.f.vec.data += self.a.harmonic_extension_trans * self.f.vec
            self.u.vec.data = solvea * self.f.vec;
            if self.static_condensation:
                self.u.vec.data += self.a.inner_solve * self.f.vec
                self.u.vec.data += self.a.harmonic_extension * self.u.vec
                
            # global last_num_its
            last_num_its = solvea.GetSteps()
            print("nze: " + str(self.a.mat.AsVector().size))
            print("number of iterations: " + str(last_num_its))
            
            
            coef_error_sqr = (self.u.components[0] - self.problemdata["Solution"])*(self.u.components[0] - self.problemdata["Solution"])
            l2diff = sqrt(IntegrateOnInterface(self.lsetmeshadap.lset_p1,self.mesh,coef_error_sqr,order=2*self.order+2))
            print("l2diff = {}".format(l2diff))
            self.mesh.UnsetDeformation()
        return l2diff

def PrintHDGTimers():            
    ### print hdg-intergrator timers
    hdgtimers = [a for a in Timers() if "HDG" in a["name"]]
    hdgtimers = sorted(hdgtimers, key=lambda k: k["time"], reverse=True)
    for timer in hdgtimers:
        print("{:<45}: {:6} cts, {:.8f} s, {:.6e} s(avg.)"
              .format(timer["name"],timer["counts"],timer["time"],timer["time"]/timer["counts"]))

import argparse        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Solve Laplace-Betrami example problem.')
    parser.add_argument('--reflvls', metavar='N', type=int, default=2, help='number of refinement levels (>=1)')
    parser.add_argument('--vtkout', dest='vtkout', action='store_true', help='enable VTKOutput')
    parser.add_argument('--no-vtkout', dest='vtkout', action='store_false', help='disable VTKOutput')
    parser.set_defaults(vtkout=False)
    args = parser.parse_args()
    options = vars(args)
    print("call arguments: ", options)

    problemdata = Make3DProblem()
    discretization = Discretization(problemdata)

    orders = [1,2,3,4,5]
    l2diffresults = [] 
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
                    l2diff = discretization.SolveProblem()
                    l2diffresults[i].append(l2diff)
                    result = True
                except Exception as e:
                    print ("Exception = ".format(e))
#                    l2diffresults[i].append(None)
                    result = False
            
            print(l2diffresults)
            print(list(map(list, zip(*l2diffresults))))
        RefineAtLevelSet(gf=discretization.lsetmeshadap.lset_p1)
    print(list(map(list, zip(*l2diffresults))))
    
    #Draw(discretization.u.components[0],discretization.mesh,"u",draw_surf=False)

    if (options['vtkout']):
        vtk = VTKOutput(ma=discretization.mesh,coefs=[discretization.lsetmeshadap.lset_p1,
                                       discretization.lsetmeshadap.deform,u.components[0]],
                        names=["lsetp1","deform","u"],filename="vtkout_",subdivision=1)
        vtk.Do()

    
#observations:
# 1. for p=1 it seems that lsetcurving p+1=2 is necessary to achieve optimal order convergence (not true for p>1)
# 2. l2order=p+1, facetorder=p gives stable result (optimal order, i.e. p+1 in the l2-norm) - however p+1 together with highest_order_dc is not stable yet (?!)

    
#[[0.8096913452313329, 0.3084025787817622, 0.11076925787924134, 0.05982871422924973], [0.11632850331154547, None, None, None], [0.03555471657060665, None, 0.00017200706247022353, None], [0.009986051916781962, 0.00018175774467425786, 7.737886013822911e-06, 3.3478037062295224e-07]]
#[[0.8096913452313312, 0.30840257878176336, 0.11076925787924342, 0.059828714229247626], [0.11632850331154539, 0.011519875896363744, 0.0013504456177305935, None], [0.03555471657059834, 0.002877418615964477, 0.00017200706247589748, 1.156914966874773e-05], [0.009986051916787501, None, 7.737886032776415e-06, None], [0.0025838732528165125, 2.249281684957255e-05, None, None]]
#[[0.8096913452313316, 0.30840257878176175, 0.11076925787923653, 0.059828714229247425], [0.11632850331153247, 0.011519875896377813, 0.001350445617732032, 0.0001810449912049936], [0.03555471657059834, 0.0028774186159388226, 0.00017200706248403185, 1.1569149801119469e-05], [0.00998605191678209, 0.00018175774468649536, 7.737886033220925e-06, 3.3478017691898923e-07]]
