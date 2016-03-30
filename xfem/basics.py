from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *

def XStdFESpace(mesh,levelset=None,flags={},basetype=None,empty=None,order=1,dirichlet=[],ref_space=None):
    """
    Creates an extended standard finite element space
    =============================================================
    an XStdFESpace is the compound fespace of a standard space 
    and an enricht space (XFESpace). The XFESpace is constrcuted
    based on the functions in the StdFESpace and a level set fct.
    =============================================================
    necessary arguments:
    - mesh : Mesh Access to mesh
    optional arguments:
    - levelset : coefficient function w.r.t. which 
                 the xfespace is created
    - flags : flags that are passed to the StdFESpace for creation
    - basetype : string defining the StdFESpace 
                 (can also be passed via flags)
    - empty : (deprecated) argument not to enrich a space 
              as a result the XFESpace will be empty but has 
              geometry information
    - order : order of StdFESpace (and hence XFESpace)
    - dirichlet : Dirichlet boundary information
    - ref_space : Subtriangulation argument 
                  choose "0" here if in doubt 

    """
    if (basetype!=None):
        flags["type_std"] = basetype
    else:
        if not ("type_std" in flags):
            flags["type_std"] = "h1ho"
    if (ref_space!=None):
        flags["ref_space"] = ref_space
    else:
        if not ("ref_space" in flags):
            flags["ref_space"] = 0
    flags["order"] = order
    if (empty!=None):
        flags["empty"] = empty
    # print("flags = {}".format(flags))
    fes = CastToXStdFESpace(FESpace("xstdfespace",mesh=mesh, flags=flags, order=order, dirichlet=dirichlet))
    if (levelset != None):
        fes.XFESpace.SetLevelSet(levelset)
        fes.Update()
        print("XStdFESpace created with {} unknowns ({} normal + {} enriched)".format(fes.ndof,fes.StdFESpace.ndof,fes.XFESpace.ndof))
    else:
        print("no levelset given, update postponed")
    return fes
        

def XFESpace(basefes,levelset=None,flags={},empty=None,ref_space=None):
    """
    Creates an extended finite element space
    =============================================================
    The XFESpace is constrcuted based on the functions in a 
    referenced base FESpace and a level set fct.
    =============================================================
    necessary arguments:
    - basefes : reference FESpace
    optional arguments:
    - levelset : coefficient function w.r.t. which 
                 the xfespace is created
    - flags : flags for the xfespace (e.g. "empty", "ref_space")
    - empty : (deprecated) argument not to enrich a space; 
              as a result the XFESpace will be empty but has 
              geometry information
    - ref_space : Subtriangulation argument 
                  choose "0" here if in doubt 

    """
    if (ref_space!=None):
        flags["ref_space"] = ref_space
    else:
        if not ("ref_space" in flags):
            flags["ref_space"] = 0
    if (empty!=None):
        flags["empty"] = empty
    fes = CastToXFESpace(FESpace("xfespace",mesh=basefes.mesh, flags=flags, order=basefes.globalorder))
    fes.SetBaseFESpace(basefes)
    if (levelset != None):
        fes.SetLevelSet(levelset)
        basefes.Update()
        fes.Update()
        print("XFESpace created with {} unknowns ({} unknowns in base FESpace)".format(fes.ndof,basefes.ndof))
    else:
        print("no levelset given, update postponed")
    return fes
        

def TwoDomainMassIntegrator (coefneg,coefpos):
    return BFI("xmass", coef=[coefneg,coefpos])

def TwoDomainSourceIntegrator (coefneg,coefpos):
    return LFI("xsource", coef=[coefneg,coefpos])



negative_domain = dict()
negative_domain["negdomain"] = True

positive_domain = dict()
positive_domain["posdomain"] = True

volume_domains = dict()
volume_domains["negdomain"] = True
volume_domains["posdomain"] = True

interface_domain = dict()
interface_domain["interface"] = True

interface_and_volume_domains = dict()
interface_and_volume_domains["negdomain"] = True
interface_and_volume_domains["posdomain"] = True
interface_and_volume_domains["interface"] = True

def IntegrateOnInterface(lset,mesh,coef,order=5,subdivlvl=0):
    return IntegrateX(lset,mesh,cf_interface=coef,order=order,
                      subdivlvl=subdivlvl,domains=interface_domain)["interface"]

def IntegrateOnPosDomain(lset,mesh,coef,order=5,subdivlvl=0):
    return IntegrateX(lset,mesh,cf_pos=coef,order=order,
                      subdivlvl=subdivlvl,domains=positive_domain)["posdomain"]

def IntegrateOnNegDomain(lset,mesh,coef,order=5,subdivlvl=0):
    return IntegrateX(lset,mesh,cf_neg=coef,order=order,
                      subdivlvl=subdivlvl,domains=negative_domain)["negdomain"]

def IntegrateOnWholeDomain(lset,mesh,cf_neg=None,cf_pos=None,coef=None,order=5,subdivlvl=0):
    if ((cf_neg == None) and (coef != None)):
        cf_neg = coef
    if ((cf_pos == None) and (coef != None)):
        cf_pos = coef
    ints = IntegrateX(lset,mesh,cf_neg=cf_neg,cf_pos=cf_pos, 
                      order=order,subdivlvl=subdivlvl,domains=volume_domains)
    return ints["negdomain"] + ints["posdomain"]


