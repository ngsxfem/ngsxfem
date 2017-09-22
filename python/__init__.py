from ngsolve.comp import *
from ngsolve.fem import *
from ngsolve.utils import L2
from xfem.ngsxfem_py import *
from xfem.ngsxfem_xfem_py import *
from xfem.ngsxfem_utils_py import *
from xfem.ngsxfem_cutint_py import *

def extend(func):
    if func.derivname == "extend":
        return func.Deriv()
    add = func.Operator("extend")
    if add:
        return add        
    raise Exception("cannot form extend")

def pos(func):
    if func.derivname == "pos":
        return func.Deriv()
    add = func.Operator("pos")
    if add:
        return add        
    raise Exception("cannot form pos")

def neg(func):
    if func.derivname == "neg":
        return func.Deriv()
    add = func.Operator("neg")
    if add:
        return add        
    raise Exception("cannot form neg")

def extend_grad(func):
    if func.derivname == "extendgrad":
        return func.Deriv()
    add = func.Operator("extendgrad")
    if add:
        return add        
    raise Exception("cannot form extend_grad")

def pos_grad(func):
    if func.derivname == "posgrad":
        return func.Deriv()
    add = func.Operator("posgrad")
    if add:
        return add        
    raise Exception("cannot form pos_grad")

def neg_grad(func):
    if func.derivname == "neggrad":
        return func.Deriv()
    add = func.Operator("neggrad")
    if add:
        return add        
    raise Exception("cannot form neg_grad")

SymbolicBFI_old = SymbolicBFI
def SymbolicBFI(levelset_domain=None, *args, **kwargs):
    if levelset_domain != None and type(levelset_domain)==dict:
        if not "force_intorder" in levelset_domain:
            levelset_domain["force_intorder"] = -1
        if not "subdivlvl" in levelset_domain:
            levelset_domain["subdivlvl"] = 0
        if not "levelset" in levelset_domain:
            print("Please provide a level set function")
        if not "domain_type" in levelset_domain:
            print("Please provide a domain type (NEG,POS or IF)")
        # print("SymbolicBFI-Wrapper: SymbolicCutBFI called")
        return SymbolicCutBFI(lset=levelset_domain["levelset"],
                              domain_type=levelset_domain["domain_type"],
                              force_intorder=levelset_domain["force_intorder"],
                              subdivlvl=levelset_domain["subdivlvl"],
                              *args, **kwargs)
    else:
        # print("SymbolicBFI-Wrapper: original SymbolicBFI called")
        if (levelset_domain == None):
            return SymbolicBFI_old(*args,**kwargs)
        else:
            return SymbolicBFI_old(levelset_domain,*args,**kwargs)

SymbolicLFI_old = SymbolicLFI
def SymbolicLFI(levelset_domain=None, *args, **kwargs):
    if levelset_domain != None and type(levelset_domain)==dict:
        if not "force_intorder" in levelset_domain:
            levelset_domain["force_intorder"] = -1
        if not "subdivlvl" in levelset_domain:
            levelset_domain["subdivlvl"] = 0
        if not "levelset" in levelset_domain:
            print("Please provide a level set function")
        if not "domain_type" in levelset_domain:
            print("Please provide a domain type (NEG,POS or IF)")
        # print("SymbolicLFI-Wrapper: SymbolicCutLFI called")
        return SymbolicCutLFI(lset=levelset_domain["levelset"],
                              domain_type=levelset_domain["domain_type"],
                              force_intorder=levelset_domain["force_intorder"],
                              subdivlvl=levelset_domain["subdivlvl"],
                              *args, **kwargs)
    else:
        # print("SymbolicLFI-Wrapper: original SymbolicLFI called")
        if (levelset_domain == None):
            return SymbolicLFI_old(*args,**kwargs)
        else:
            return SymbolicLFI_old(levelset_domain,*args,**kwargs)

def Integrate_X_special_args(levelset_domain={}, cf=None, mesh=None, VOL_or_BND=VOL, order=5, region_wise=False, element_wise = False, heapsize=1000000):
    if not "force_intorder" in levelset_domain or levelset_domain["force_intorder"] == -1:
        levelset_domain["force_intorder"] = -1
    else:
        order = levelset_domain["force_intorder"]
        
    if not "subdivlvl" in levelset_domain:
        levelset_domain["subdivlvl"] = 0
    if not "levelset" in levelset_domain:
        print("Please provide a level set function")
    if not "domain_type" in levelset_domain:
        print("Please provide a domain type (NEG,POS or IF)")

    return IntegrateX(lset=levelset_domain["levelset"],
                      mesh=mesh, cf=cf,
                      order=order,
                      domain_type=levelset_domain["domain_type"],
                      subdivlvl=levelset_domain["subdivlvl"],
                      heapsize=heapsize)
    

##### THIS IS ANOTHER WRAPPER (original IntegrateX-interface is pretty ugly...) TODO        
Integrate_old = Integrate
def Integrate(levelset_domain=None, *args, **kwargs):
    if levelset_domain != None and type(levelset_domain)==dict:
        # print("Integrate-Wrapper: IntegrateX called")
        return Integrate_X_special_args(levelset_domain, *args, **kwargs)
    else:
        # print("Integrate-Wrapper: original Integrate called")
        if (levelset_domain == None):
            return Integrate_old(*args,**kwargs)
        else:
            newargs = [levelset_domain]
            for q in args:
                newargs.append(q)
            return Integrate_old(*newargs,**kwargs)

def IndicatorCF(mesh, ba, facets = False):
    if facets:
        ret = GridFunction(FESpace("facet",mesh,order=0))
        for i in range(len(ba)):
            if ba[i]:
                ret.vec[i] = 1.0
            else:
                ret.vec[i] = 0.0
        return ret
    else:
        return BitArrayCF(BitArray(ba))

def CutRatioGF(cutinfo):
    ret = GridFunction(L2(cutinfo.Mesh(),order=0))
    ret.vec.data = cutinfo.GetCutRatios(VOL)
    return ret

def kappa(mesh,lset_approx, subdivlvl=0):
    print("kappa-function is deprecated - use CutRatioGF instead")
    kappa1 = GridFunction(L2(mesh,order=0))
    lset_neg = { "levelset" : lset_approx, "domain_type" : NEG, "subdivlvl" : subdivlvl}
    kappa_f = LinearForm(kappa1.space)
    kappa_f += SymbolicLFI(levelset_domain = lset_neg, form = kappa1.space.TestFunction() )
    kappa_f.Assemble();
    kappa1.space.SolveM(CoefficientFunction(1.0),kappa_f.vec)
    kappa1.vec.data = kappa_f.vec
    kappa2 = 1.0 - kappa1
    return (kappa1,kappa2)

def IsCut(mesh,lset_approx, subdivlvl=0):
    def cut(kappa):
        if (kappa > 1e-16 and kappa < 1.0-1e-16):
            return 1.0
        else:
            return 0.0
    kappa1, dummycf = kappa(mesh,lset_approx,subdivlvl)
    for i in range(len(kappa1.vec)):
        kappa1.vec[i] = cut(kappa1.vec[i])
    # from numpy import vectorize
    # cut = vectorize(cut)
    # kappa1.vec.FV().NumPy()[:] = cut(kappa1.vec.FV().NumPy())
    return kappa1

