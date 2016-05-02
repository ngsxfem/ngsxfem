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
        
def TwoDomainLaplaceIntegrator (coefneg,coefpos):
    return BFI("xlaplace", coef=[CoefficientFunction(coefneg),CoefficientFunction(coefpos)])

def TwoDomainMassIntegrator (coefneg,coefpos):
    return BFI("xmass", coef=[CoefficientFunction(coefneg),CoefficientFunction(coefpos)])

def TwoDomainSourceIntegrator (coefneg,coefpos):
    return LFI("xsource", coef=[CoefficientFunction(coefneg),CoefficientFunction(coefpos)])

def GhostPenaltyIntegrator (coefneg=1.0,coefpos=1.0,stab_param=1.0):
    return BFI("lo_ghostpenalty", coef=[CoefficientFunction(coefneg),CoefficientFunction(coefpos),CoefficientFunction(stab_param)])



def XNitscheIntegrators (diffusion, henryweights=[1.0,1.0], weighting="hansbo", stab_param=10.0, minstab=False,fluxjump=None,jump=None):
    """
    Returns Integrators for (scalar) Nitsche-XFEM discretiz.,
    i.e. weak enforcement of interface conditions in XFEM-setting
    interface conditions are:
    [u] = g 
    and [a du/dn] = h with a the piecewise constant diffusion
    =============================================================
    Integrators for (scalar) Nitsche-XFEM discretizations    
    returns a 3-tuple of Integrators: 
      1. BFI: NitscheXFEM - left hand side
      2. LFI: NitscheXFEM - r.h.s. w.r.t. jump g in [u] = g
      3. LFI: NitscheXFEM - r.h.s. w.r.t. jump h in [a du/dn] = h
    =============================================================
    necessary arguments:
    - diffusion : 
        list of diffusion constants
    - henryweights : 
        list of henry weights (default = [1,1])
    - weighting : 
        weighting of the averaging in  Nitsche form.
        {u} = k_1 u_1 + k_2 * u_2
        possibilities are (k_2 = 1.0 - k_1):
         - "hansbo" (default): 
            volume weighted average:
              k_1 = |T_1| / |T|
         - "halfhalf" / "naive": 
            equal weights:
              k_1 = k_2 = 0.5
         - "heaviside": 
              k_1 = 1.0 if |T_1| / |T| > 0.5
              k_1 = 0.0 else
        default is "hansbo"
    - stab_param : 
        penalty parameter in penalty of Nitsche form.
        default is 10.0.
    - minstab : 
        Formulation with lifting stabilization. 
        No stab_param needed in this formulation (value ignored).
        default = False
    - jump : 
        coefficient g in interface condition [u] = h
        default is None, i.e. [u] = 0
    - fluxjump : 
        coefficient h in interface condition [a du/dn] = h
        default is None, i.e. [a du/dn] = 0
    """
    
    rhs1 = None
    rhs2 = None
    
    if minstab:
        coefs=[diffusion[0], diffusion[1], henryweights[0], henryweights[1]]
        if fluxjump != None:
            coefs_fluxjump=[CoefficientFunction(henryweights[0]),
                            CoefficientFunction(henryweights[1]),
                            CoefficientFunction(fluxjump)]
        if jump != None:
            coefs_jump=[CoefficientFunction(diffusion[0]),
                        CoefficientFunction(diffusion[1]),
                        CoefficientFunction(henryweights[0]),
                        CoefficientFunction(henryweights[1]),
                        CoefficientFunction(jump)]

        if weighting.lower() == "hansbo":
            lhs = BFI("xnitsche_minstab_hansbo", coef=coefs)
            if jump != None:
                rhs1 = LFI("xnitscherhsjump_hansbo", coef=coefs_jump)
            if fluxjump != None:
                rhs2 = LFI("xnitscherhsfluxjump_hansbo", coef=coefs_fluxjump)
        elif weighting.lower() == "halfhalf" or weighting.lower() == "naive":
            lhs = BFI("xnitsche_minstab_halfhalf", coef=coefs)
            if jump != None:
                rhs1 = LFI("xnitscherhsjump_halfhalf", coef=coefs_jump)
            if fluxjump != None:
                rhs2 = LFI("xnitscherhsfluxjump_halfhalf", coef=coefs_fluxjump)
        elif weighting.lower() == "heaviside":
            lhs = BFI("xnitsche_minstab_heaviside", coef=coefs)
            if jump != None:
                rhs1 = LFI("xnitscherhsjump_heaviside", coef=coefs_jump)
            if fluxjump != None:
                rhs2 = LFI("xnitscherhsfluxjump_heaviside", coef=coefs_fluxjump)
        else:
            print ("unknown weighting for nitsche integrators!")
    else:
        coefs=[diffusion[0], diffusion[1], henryweights[0], henryweights[1], stab_param]
        if fluxjump != None:
            coefs_fluxjump=[CoefficientFunction(henryweights[0]),
                            CoefficientFunction(henryweights[1]),
                            CoefficientFunction(fluxjump)]
        if jump != None:
            coefs_jump=[CoefficientFunction(diffusion[0]),
                        CoefficientFunction(diffusion[1]),
                        CoefficientFunction(henryweights[0]),
                        CoefficientFunction(henryweights[1]),
                        CoefficientFunction(jump),
                        CoefficientFunction(stab_param)]
            
        if weighting.lower() == "hansbo":
            lhs = BFI("xnitsche_hansbo", coef=coefs)
            if jump != None:
                rhs1 = LFI("xnitscherhsjump_hansbo", coef=coefs_jump)
            if fluxjump != None:
                rhs2 = LFI("xnitscherhsfluxjump_hansbo", coef=coefs_fluxjump)
        elif weighting.lower() == "halfhalf" or weighting.lower() == "naive":
            lhs = BFI("xnitsche_halfhalf", coef=coefs)
            if jump != None:
                rhs1 = LFI("xnitscherhsjump_halfhalf", coef=coefs_jump)
            if fluxjump != None:
                rhs2 = LFI("xnitscherhsfluxjump_halfhalf", coef=coefs_fluxjump)
        elif weighting.lower() == "heaviside":
            lhs = BFI("xnitsche_heaviside", coef=coefs)
            if jump != None:
                rhs1 = LFI("xnitscherhsjump_heaviside", coef=coefs_jump)
            if fluxjump != None:
                rhs2 = LFI("xnitscherhsfluxjump_heaviside", coef=coefs_fluxjump)
        else:
            print ("unknown weighting for nitsche integrators!")
        
    return (lhs,rhs1,rhs2)


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

def IntegrateOnInterface(lset,mesh,coef,order=5,subdivlvl=0,heapsize=1000000):
    return IntegrateX(lset,mesh,cf_interface=coef,order=order,
                      subdivlvl=subdivlvl,domains=interface_domain,heapsize=heapsize)["interface"]

def IntegrateOnPosDomain(lset,mesh,coef,order=5,subdivlvl=0,heapsize=1000000):
    return IntegrateX(lset,mesh,cf_pos=coef,order=order,
                      subdivlvl=subdivlvl,domains=positive_domain,heapsize=heapsize)["posdomain"]

def IntegrateOnNegDomain(lset,mesh,coef,order=5,subdivlvl=0,heapsize=1000000):
    return IntegrateX(lset,mesh,cf_neg=coef,order=order,
                      subdivlvl=subdivlvl,domains=negative_domain,heapsize=heapsize)["negdomain"]

def IntegrateOnWholeDomain(lset,mesh,cf_neg=None,cf_pos=None,coef=None,order=5,subdivlvl=0,heapsize=1000000):
    if ((cf_neg == None) and (coef != None)):
        cf_neg = coef
    if ((cf_pos == None) and (coef != None)):
        cf_pos = coef
    ints = IntegrateX(lset,mesh,cf_neg=cf_neg,cf_pos=cf_pos, 
                      order=order,subdivlvl=subdivlvl,domains=volume_domains,heapsize=heapsize)
    return ints["negdomain"] + ints["posdomain"]


