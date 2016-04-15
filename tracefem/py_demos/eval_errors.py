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

def EvalErrorsOnInterface(mesh,lsetmeshadap,discrete_solution,order,problemdata,resultdict):
    levelset_p1 = lsetmeshadap.lset_p1

    mesh.SetDeformation(lsetmeshadap.deform)

    if problemdata["Solution"] != None:

        if (problemdata["Reaction"] == None or problemdata["Reaction"] == 0.0):
            surface_meas = IntegrateOnInterface(lsetmeshadap.lset_p1,mesh,CoefficientFunction(1.0),order=order,heapsize=10000000)
            total_u = IntegrateOnInterface(lsetmeshadap.lset_p1,mesh,CoefficientFunction(discrete_solution)-problemdata["Solution"],order=order,heapsize=10000000)
            correction_constant = total_u / surface_meas
            coef_error_sqr = (discrete_solution - correction_constant - problemdata["Solution"])*(discrete_solution - correction_constant - problemdata["Solution"])
        else:
            coef_error_sqr = (discrete_solution - problemdata["Solution"])*(discrete_solution - problemdata["Solution"])

        l2diff = sqrt(IntegrateOnInterface(levelset_p1,mesh,coef_error_sqr,order=2*order+2,heapsize=10000000))
        print("l2diff = {}".format(l2diff))
        resultdict["l2err"] = l2diff
        
    nhelp = levelset_p1.Deriv()
    n = 1.0/sqrt(nhelp*nhelp) * nhelp
    un = discrete_solution.Deriv()*n

    coef_gradnormal_error = un*un
    h1normdiff = sqrt(IntegrateOnInterface(levelset_p1,mesh,coef_gradnormal_error,order=2*order,heapsize=10000000))
    print("h1normerr = {}".format(h1normdiff))
    resultdict["h1normerr"] = h1normdiff
        
    if problemdata["GradSolution"] != None:
        tanggrad = discrete_solution.Deriv() - un * n
        coef_gradtang_error = (tanggrad - problemdata["GradSolution"])*(tanggrad - problemdata["GradSolution"])
        h1tangdiff = sqrt(IntegrateOnInterface(levelset_p1,mesh,coef_gradtang_error,order=2*order,heapsize=10000000))
        print("h1tangerr = {}".format(h1tangdiff))
        resultdict["h1tangerr"] = h1tangdiff
        
    maxdistlset = lsetmeshadap.CalcMaxDistance(problemdata["Levelset"],heapsize=10000000);
    print("maxdist = {}".format(maxdistlset))
    resultdict["maxdist"] = maxdistlset
        
    mesh.UnsetDeformation()

