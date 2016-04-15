# ngsolve stuff
from ngsolve import *
# xfem and trace fem stuff
import libngsxfem_py.xfem as xfem                                 
# For Integration on Interface
from xfem.basics import *
# For LevelSetAdaptationMachinery
from xfem.lsetcurv import *
# For TraceFEM-Integrators (convenience)
from xfem.tracefem import *


def constraint_pcg(mat, pre, rhs, constraint, e, prec = 1e-8, maxits = 100):
    """preconditioned conjugate gradient method"""

    u = rhs.CreateVector()
    d = rhs.CreateVector()
    w = rhs.CreateVector()
    s = rhs.CreateVector()

    u *= 0
    d.data = rhs - mat * u
    w.data = pre * d
    s.data = w
    wdn = InnerProduct (w,d)
    print ("|u0| = ", wdn)
    
    for it in range(maxits):
        w.data = mat * s
        wd = wdn
        Ass = InnerProduct (s, w)
        alpha = wd / Ass
        u.data += alpha * s
        d.data += (-alpha) * w

        w.data = pre * d
        wdn = InnerProduct (w, d)
        beta = wdn / wd

        s *= beta
        s.data += w

        err = sqrt(wd)
        if err < prec:
            break
        print ("\rit = ", it, " err = ", err, end="")

        corr = InnerProduct(constraint,u)
        u.data -= corr * e

    print("")
    return (u,it)






class TraceFEMDiscretization(object):
    """
    Tracefem(continuous)Discretization: definition of FESpaces, bi/linear forms, etc..
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
            
        self.lsetmeshadap = LevelSetMeshAdaptation(self.mesh, order=self.order+1, threshold=1000, discontinuous_qn=True,heapsize=10000000)

        symmetric = True
        if (problemdata["Convection"] != None):
            symmetric = False
            

        ### Setting up discrete variational problem
        self.Vh = H1(self.mesh, order=self.order, dirichlet=[])
        self.Vh_tr = TraceFESpace(self.mesh, self.Vh, problemdata["Levelset"])

        self.a = BilinearForm(self.Vh_tr, symmetric = symmetric, flags = {"eliminate_internal" : self.static_condensation})
        if (problemdata["Reaction"] != None):
            self.a += TraceMass(problemdata["Reaction"])
        if (problemdata["Diffusion"] != None):
            self.a += TraceLaplaceBeltrami(problemdata["Diffusion"])
        if (problemdata["VolumeStabilization"]!=None):
            self.a += NormalLaplaceStabilization(problemdata["VolumeStabilization"],
                                                 # problemdata["GradLevelset"])
                                                 self.lsetmeshadap.lset_p1.Deriv())
        if (problemdata["Convection"] != None):
            self.a += TraceConvection(problemdata["Convection"])

        if (self.problemdata["Iterative"]):
            self.c = Preconditioner(self.a, type="local")
        else:
            self.c = Preconditioner(self.a, type="direct")
            
        self.u = GridFunction(self.Vh_tr)
  
    def UpdateSpace(self):
        self.Vh.Update()
        self.Vh_tr.Update()

    def VolumeSolution(self):
        return self.u
                    
    def SolveProblem(self,firstcall=False, results = None):
        if results == None:
            results = {}
        # Calculation of the deformation:
        deformation = self.lsetmeshadap.CalcDeformation(self.problemdata["Levelset"])
        # Applying the mesh deformation
        self.mesh.SetDeformation(deformation)

        if (not firstcall):
            self.UpdateSpace()
        
        results["total_ndofs"] = self.Vh_tr.ndof
        global_ndofs = self.Vh_tr.ndof
        for i in range(self.Vh_tr.ndof):
        #     print (str(i) + " : " + str(Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF))
            if (self.Vh_tr.CouplingType(i) == COUPLING_TYPE.LOCAL_DOF):
                # self.Vh_tr.FreeDofs()[i] = 0
                global_ndofs = global_ndofs - 1
        results["global_ndofs"] = global_ndofs
        
        self.u.Update()

        self.f = LinearForm(self.Vh_tr)
        if ("SourceMeanValue" in self.problemdata and self.problemdata["SourceMeanValue"] != None):
            surface_meas = IntegrateOnInterface(self.lsetmeshadap.lset_p1,self.mesh,CoefficientFunction(1.0),order=self.order,heapsize=10000000)
            total_f = IntegrateOnInterface(self.lsetmeshadap.lset_p1,self.mesh,self.problemdata["Source"],order=self.order,heapsize=10000000)
            f_h = self.problemdata["Source"] + ( self.problemdata["SourceMeanValue"] - total_f/surface_meas )
            self.f += TraceSource(f_h)
            print("correction constant in r.h.s. (to match SourceMeanValue) = {}".format(total_f))
        else:
            self.f += TraceSource(self.problemdata["Source"])

        constraint = LinearForm(self.Vh_tr)
        constraint += TraceSource(CoefficientFunction(1.0))
        constraint.Assemble(heapsize=10000000);

        surface_meas = IntegrateOnInterface(self.lsetmeshadap.lset_p1,self.mesh,CoefficientFunction(1.0),order=self.order,heapsize=10000000)

        m = BilinearForm(self.Vh_tr)
        m += TraceMass(CoefficientFunction(surface_meas))
        m += NormalLaplaceStabilization(specialcf.mesh_size,self.lsetmeshadap.lset_p1.Deriv())

        m.Assemble(heapsize=10000000);

        e = self.u.vec.CreateVector()
        e.data = m.mat.Inverse() * constraint.vec
        
        print (e.data)
        input("q")
        print ("Asf = ", InnerProduct(e,constraint.vec))

        self.a.Assemble(heapsize=10000000);
        self.f.Assemble(heapsize=10000000);
        self.c.Update();

        self.u.vec[:] = 0
        
        # pcg(self.a.mat, self.c.pre, self.f.vec, maxits = 100)
        # solvea = CGSolver( mat=self.a.mat, pre=self.c.mat, complex=False, printrates=False, precision=1e-8, maxsteps=2000000)
        # u.vec.data = ainv * f.vec
        
        if self.static_condensation:
            self.f.vec.data += self.a.harmonic_extension_trans * self.f.vec
        self.u.vec.data, last_num_its = constraint_pcg(self.a.mat, self.c.mat, self.f.vec, constraint.vec, e, prec=abs(total_f), maxits = 10000)
        # self.u.vec.data = solvea * self.f.vec;
        if self.static_condensation:
            self.u.vec.data += self.a.inner_solve * self.f.vec
            self.u.vec.data += self.a.harmonic_extension * self.u.vec
            
        print("nze: " + str(self.a.mat.AsVector().size))
        results["nze"] = self.a.mat.AsVector().size
        
        # last_num_its = solvea.GetSteps()
        print("number of iterations: " + str(last_num_its))
        results["numits"] = last_num_its
        
        self.mesh.UnsetDeformation()
        return results

