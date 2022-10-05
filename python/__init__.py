"""
(ngs)xfem
=========

A module for unfitted finite element discretizations in NGSolve

Submodules:
xfem.cutmg ... MultiGrid for CutFEM
xfem.lsetcurving ... isoparametric unfitted FEM
xfem.lset_spacetime ... isoparametric unfitted space-time FEM
xfem.mlset ... multiple level sets
xfem.utils ... some example level set geometries
"""



from ngsolve import (L2, VOL, BitArray, CoefficientFunction, FESpace,
                     GridFunction, H1, IfPos, LinearForm, Parameter, dx)
from ngsolve.comp import Integrate as ngsolve_Integrate
from ngsolve.comp import ProxyFunction
from ngsolve.comp import SymbolicBFI as ngsolve_SymbolicBFI
from ngsolve.comp import SymbolicLFI as ngsolve_SymbolicLFI
from xfem.ngsxfem_py import *
from xfem.ngs_check import check_if_ngsolve_newer_than, __ngsolve_required__

check_if_ngsolve_newer_than(__ngsolve_required__)

def HAS(domain_type):
    """
For a given domain_type return the combined domain type that 
includes all elements that have a part in the domain type.
    """
    if domain_type == NEG:
        return HASNEG
    elif domain_type == POS:
        return HASPOS
    elif domain_type == IF:
        return IF
    else:
        raise Exception("invalid domain type")

def extend(func):
    """
Evaluates the XFiniteElement-function independent of the level set domains.

Note:
This will lead to the same behavior as the function that the XFiniteElement-function is based
on.
    """
    if func.derivname == "extend":
        return func.Deriv()
    add = func.Operator("extend")
    if add:
        return add
    raise Exception("cannot form extend")

def pos(func):
    """
Evaluates an XFiniteElement-function assuming a positive level set domain.

Note:
This can lead to non-zero values also in domains where the level set function is non-positive.
    """
    if func.derivname == "pos":
        return func.Deriv()
    add = func.Operator("pos")
    if add:
        return add
    raise Exception("cannot form pos")

def neg(func):
    """
Evaluates an XFiniteElement-function assuming a negative level set domain.

Note:
This can lead to non-zero values also in domains where the level set function is non-negative.
    """
    if func.derivname == "neg":
        return func.Deriv()
    add = func.Operator("neg")
    if add:
        return add
    raise Exception("cannot form neg")

def extend_grad(func):
    """
Evaluates the gradient of an XFiniteElement-function independent of the level set domains.

Note:
This will lead to the same behavior as the function that the XFiniteElement-function is based on.
    """
    if func.derivname == "extendgrad":
        return func.Deriv()
    add = func.Operator("extendgrad")
    if add:
        return add
    raise Exception("cannot form extend_grad")

def pos_grad(func):
    """
Evaluates the gradient of an XFiniteElement-function assuming a positive level set domain.

Note:
This can lead to non-zero values also in domains where the level set function is non-positive.
    """
    if func.derivname == "posgrad":
        return func.Deriv()
    add = func.Operator("posgrad")
    if add:
        return add
    raise Exception("cannot form pos_grad")

def neg_grad(func):
    """
Evaluates the gradient of an XFiniteElement-function assuming a negative level set domain.

Note:
This can lead to non-zero values also in domains where the level set function is non-negative.
    """
    if func.derivname == "neggrad":
        return func.Deriv()
    add = func.Operator("neggrad")
    if add:
        return add
    raise Exception("cannot form neg_grad")

def dtref(func):
    """
Evaluates the time derivative (w.r.t. the reference time interval) of a Space-Time function.
    """
    add = func.Operator("dt")
    if add:
        return add
    raise Exception("cannot form dt")

def dt(func):
    """
Deprecated: use "dtref" instead
    """
    print("WARNING: dt is deprecated. Use \"dtref\" instead. \n         Note that the operator acts w.r.t. the reference time intervals.")
    return dtref(func)



def SymbolicBFIWrapper(levelset_domain=None, *args, **kwargs):
    """
Wrapper around SymbolicBFI to allow for integrators on level set domains (see also
SymbolicCutBFI). The dictionary contains the level set function (CoefficientFunciton or
GridFunction) and the domain-type (NEG/POS/IF). If the dictionary is not provided, the standard
SymbolicBFI function from NGSolve will be called.

Parameters

levelset_domain : dictionary
  entries:
  * "levelset": 
    singe level set : ngsolve.CoefficientFunction
      CoefficientFunction that describes the geometry. In the best case lset is a GridFunction of an
      FESpace with scalar continuous piecewise (multi-) linear basis functions.
    multiple level sets: tuple(ngsolve.GridFunction)
      Tuple of GridFunctions that describe the geometry.
  * "domain_type" :
    single level set: {NEG,POS,IF} (ENUM) 
      Integration on the domain where either:
      * the level set function is negative (NEG)
      * the level set function is positive (POS)
      * the level set function is zero     (IF )
    multiple level sets: {tuple({ENUM}), list(tuple(ENUM)), DomainTypeArray}
      Integration on the domains specified
  * "subdivlvl" : int
    On simplex meshes a subtriangulation is created on which the level set function lset is
    interpolated piecewise linearly. Based on this approximation, the integration rule is
    constructed. Note: this argument only works on simplices.
  * "order" : int
    (default: entry does not exist or value -1)
    overwrites "order"-arguments in the integration
  * "quad_dir_policy" : {FIRST, OPTIMAL, FALLBACK} (ENUM)
    Integration direction policy for iterated integrals approach
    * first direction is used unless not applicable (FIRST)
    * best direction (in terms of transformation constant) is used (OPTIMAL)
    * subdivision into simplices is always used (FALLBACK)

Other Parameters :

  form : ngsolve.CoefficientFunction
    form to integrate

  VOL_or_BND : {VOL,BND}
    integrator is defined in the volume or boundary

  element_boundary : boolean
    Integration of the boundary of an element
    (not active for level set domains)

  skeleton : boolean
    Integration over element-interface

  definedon : Region
    Domain description on where the integrator is defined

  definedonelements: BitArray
    BitArray that allows integration only on elements or facets (if skeleton=True) that are marked
    True.

  deformation : GridFunction
    Specify a specific mesh deformation for a bilinear form

  order : int
    Modifies the order of the quadrature rule used. This is overruled by "order"-entry of the 
    levelset_domain dictionary, if the dictionary entry exists.

  time_order : int
    order in time that is used in the space-time integration. time_order=-1 means that no space-time
    rule will be applied. This is only relevant for space-time discretizations.
"""
    if levelset_domain != None and type(levelset_domain)==dict:
        # shallow copy is sufficient to modify "order" and "time_order" locally
        levelset_domain_local = levelset_domain.copy()
        if "order" in kwargs:
            if not "order" in levelset_domain_local:
                levelset_domain_local["order"] = kwargs["order"]
            del kwargs["order"]
        if "time_order" in kwargs:
            if not "time_order" in levelset_domain_local or levelset_domain_local["time_order"] == -1:
                levelset_domain_local["time_order"] = kwargs["time_order"]
            del kwargs["time_order"]

        # print("SymbolicBFI-Wrapper: SymbolicCutBFI called")
        return SymbolicCutBFI(levelset_domain=levelset_domain_local,*args, **kwargs)
    else:
        # print("SymbolicBFI-Wrapper: original SymbolicBFI called")
        if (levelset_domain == None):
            return ngsolve_SymbolicBFI(*args,**kwargs)
        else:
            return ngsolve_SymbolicBFI(levelset_domain,*args,**kwargs)

def SymbolicLFIWrapper(levelset_domain=None, *args, **kwargs):
    """
Wrapper around SymbolicLFI to allow for integrators on level set domains (see also
SymbolicCutLFI). The dictionary contains the level set function (CoefficientFunciton or
GridFunction) and the domain-type (NEG/POS/IF). If the dictionary is not provided, the standard
SymbolicLFI function from NGSolve will be called.

Parameters

levelset_domain : dictionary
  entries:
  * "levelset": 
    singe level set : ngsolve.CoefficientFunction
      CoefficientFunction that describes the geometry. In the best case lset is a GridFunction of an
      FESpace with scalar continuous piecewise (multi-) linear basis functions.
    multiple level sets: tuple(ngsolve.GridFunction)
      Tuple of GridFunctions that describe the geometry.
  * "domain_type" :
    single level set: {NEG,POS,IF} (ENUM) 
      Integration on the domain where either:
      * the level set function is negative (NEG)
      * the level set function is positive (POS)
      * the level set function is zero     (IF )
    multiple level sets: {tuple({ENUM}), list(tuple(ENUM)), DomainTypeArray}
      Integration on the domains specified
  * "subdivlvl" : int
    On simplex meshes a subtriangulation is created on which the level set function lset is
    interpolated piecewise linearly. Based on this approximation, the integration rule is
    constructed. Note: this argument only works on simplices.
  * "order" : int
    (default: entry does not exist or value -1)
    overwrites "order"-arguments in the integration
  * "quad_dir_policy" : {FIRST, OPTIMAL, FALLBACK} (ENUM)
    Integration direction policy for iterated integrals approach
    * first direction is used unless not applicable (FIRST)
    * best direction (in terms of transformation constant) is used (OPTIMAL)
    * subdivision into simplices is always used (FALLBACK)

Other Parameters :

  form : ngsolve.CoefficientFunction
    form to integrate

  VOL_or_BND : {VOL,BND}
    integrator is defined in the volume or boundary

  element_boundary : boolean
    Integration of the boundary of an element
    (not active for level set domains)

  skeleton : boolean
    Integration over element-interface

  definedon : Region
    Domain description on where the integrator is defined

  definedonelements: BitArray
    BitArray that allows integration only on elements or facets (if skeleton=True) that are marked
    True.

  deformation : GridFunction
      Specify a specific mesh deformation for a linear form

  order : int
    Modifies the order of the quadrature rule used. This is overruled by "order"-entry of the 
    levelset_domain dictionary, if the dictionary entry exists.

  time_order : int
    order in time that is used in the space-time integration. time_order=-1 means that no space-time
    rule will be applied. This is only relevant for space-time discretizations. Note that
    time_order can only be active if the key "time_order" of the levelset_domain is not set (or -1)
"""
    if levelset_domain != None and type(levelset_domain)==dict:
        # shallow copy is sufficient to modify "order" and "time_order" locally
        levelset_domain_local = levelset_domain.copy()      
        if "order" in kwargs:
            if not "order" in levelset_domain_local:
                levelset_domain_local["order"] = kwargs["order"]
            del kwargs["order"]
        if "time_order" in kwargs:
            if not "time_order" in levelset_domain_local or levelset_domain_local["time_order"] == -1:
                levelset_domain_local["time_order"] = kwargs["time_order"]
            del kwargs["time_order"]
        
        return SymbolicCutLFI(levelset_domain=levelset_domain_local,*args, **kwargs)
    else:
        if (levelset_domain == None):
            return ngsolve_SymbolicLFI(*args,**kwargs)
        else:
            return ngsolve_SymbolicLFI(levelset_domain,*args,**kwargs)

def Integrate_X_special_args(levelset_domain={}, cf=None, mesh=None, VOL_or_BND=VOL, order=5, time_order=-1, region_wise=False, element_wise = False, heapsize=1000000, ip_container=None):
    """
Integrate_X_special_args should not be called directly.
See documentation of Integrate.
    """
    levelset_domain_local = levelset_domain.copy() 
    if not "order" in levelset_domain_local or levelset_domain_local["order"] == -1:
        levelset_domain_local["order"] = order
    if not "time_order" in levelset_domain_local or levelset_domain_local["time_order"] == -1:
        levelset_domain_local["time_order"] = time_order
    return IntegrateX(levelset_domain = levelset_domain_local,
                      mesh=mesh, cf=cf,
                      ip_container=ip_container,
                      element_wise=element_wise,
                      heapsize=heapsize
                      )


##### THIS IS ANOTHER WRAPPER (original IntegrateX-interface is pretty ugly...) TODO
def Integrate(levelset_domain=None, *args, **kwargs):
    """
Integrate-wrapper. If a dictionary 'levelset_domain' is provided integration will be done on the
level set part of the mesh. The dictionary contains the level set function (CoefficientFunciton or
GridFunction) and the domain-type (NEG/POS/IF). If the dictionary is not provided, the standard
Integrate function from NGSolve will be called.

Parameters

levelset_domain : dictionary
  entries:
  * "levelset": 
    singe level set : ngsolve.CoefficientFunction
      CoefficientFunction that describes the geometry. In the best case lset is a GridFunction of an
      FESpace with scalar continuous piecewise (multi-) linear basis functions.
    multiple level sets: tuple(ngsolve.GridFunction)
      Tuple of GridFunctions that describe the geometry.
  * "domain_type" :
    single level set: {NEG,POS,IF} (ENUM) 
      Integration on the domain where either:
      * the level set function is negative (NEG)
      * the level set function is positive (POS)
      * the level set function is zero     (IF )
    multiple level sets: {tuple({ENUM}), list(tuple(ENUM)), DomainTypeArray}
      Integration on the domains specified
  * "subdivlvl" : int
    On simplex meshes a subtriangulation is created on which the level set function lset is
    interpolated piecewise linearly. Based on this approximation, the integration rule is
    constructed. Note: this argument only works on simplices.
  * "order" : int
    (default: entry does not exist or value -1)
    overwrites "order"-arguments in the integration (affects only spatial integration)
  * "time_order" : int 
    defines integration order in time (for space-time integrals only)
  * "quad_dir_policy" : {FIRST, OPTIMAL, FALLBACK} (ENUM)
    Integration direction policy for iterated integrals approach
    * first direction is used unless not applicable (FIRST)
    * best direction (in terms of transformation constant) is used (OPTIMAL)
    * subdivision into simplices is always used (FALLBACK)

mesh :
  Mesh to integrate on (on some part)

cf : ngsolve.CoefficientFunction
  the integrand

order : int (default = 5)
  integration order. Can be overruled by "order"-entry of the levelset_domain dictionary.

time_order : int (default = -1)
  integration order in time (for space-time integration), default: -1 (no space-time integrals)

region_wise : bool
  (only active for non-levelset version)

element_wise : bool
  integration result is return per element

ip_container : list (or None)
  a list to store integration points (for debugging or visualization purposes only!)

heapsize : int
  heapsize for local computations.
    """
    if levelset_domain != None and type(levelset_domain)==dict:
        # print("Integrate-Wrapper: IntegrateX called")
        return Integrate_X_special_args(levelset_domain, *args, **kwargs)
    else:
        # print("Integrate-Wrapper: original Integrate called")
        if (levelset_domain == None):
            return ngsolve_Integrate(*args,**kwargs)
        else:
            newargs = [levelset_domain]
            for q in args:
                newargs.append(q)
            return ngsolve_Integrate(*newargs,**kwargs)

def IndicatorCF(mesh, ba, facets = False):
    """
Returns a CoefficientFunction that evaluates a BitArray. On elements/facets with an index i where
the BitArray evaluates to true the CoefficientFunction will evaluate as 1, otherwise as 0. Similar
functionality (only on elements) can be obtained with BitArrayCF.
    """
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
    """
Ratio between negative and full part of an element. Vector taken from CutInfo and put into a
piecewise constant GridFunction.
    """
    ret = GridFunction(L2(cutinfo.Mesh(),order=0))
    ret.vec.data = cutinfo.GetCutRatios(VOL)
    return ret

def kappa(mesh,lset_approx, subdivlvl=0):
    """
Tuple of ratios between negative/positive and full
part of an element (deprecated).
    """
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
    """
GridFunction that is 1 on cut elements, 0 otherwise (deprecated). Use CutInfo-functionality (perhaps
combined with BitArrayCF).
    """
    print("IsCut-function is deprecated - use CutInfo-functionality instead")
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

all_domain_types = [ DOMAIN_TYPE.NEG,
                     DOMAIN_TYPE.POS,
                     DOMAIN_TYPE.IF ]            

all_combined_domain_types = [ COMBINED_DOMAIN_TYPE.NO,
                              COMBINED_DOMAIN_TYPE.CDOM_NEG,
                              COMBINED_DOMAIN_TYPE.CDOM_POS,
                              COMBINED_DOMAIN_TYPE.UNCUT,
                              COMBINED_DOMAIN_TYPE.CDOM_IF,
                              COMBINED_DOMAIN_TYPE.HASNEG,
                              COMBINED_DOMAIN_TYPE.HASPOS,
                              COMBINED_DOMAIN_TYPE.ANY ]            

def SpaceTimeWeakSet(gfu_e, cf, space_fes):
    """
Ondocumented feature
    """
    gfu_e_repl = GridFunction(space_fes)
    gfu_e_repl.Set( cf )
    gfu_e.vec[:].data = gfu_e_repl.vec

ngsolveSet = GridFunction.Set
def SpaceTimeSet(self, cf, *args, **kwargs):
    """
Overrides the NGSolve version of Set in case of a space-time FESpace.
In this case the usual Set() is used on each nodal dof in time.
    """
    if (isinstance(self.space,CSpaceTimeFESpace)):
      cf = CoefficientFunction(cf)
      gfs = GridFunction(self.space.spaceFES)
      ndof_node = len(gfs.vec)
      j = 0
      for i,ti in enumerate(self.space.TimeFE_nodes()):
        if self.space.IsTimeNodeActive(i):
          ngsolveSet(gfs,fix_tref(cf,ti), *args, **kwargs)
          self.vec[j*ndof_node : (j+1)*ndof_node].data = gfs.vec[:]
          j += 1
    else:
      ngsolveSet(self,cf, *args, **kwargs)

def fix_tref(obj,time,*args,**kwargs):
    """
Takes a (possibly space-time) CoefficientFunction and fixes the temporal
variable to `time` and return this as a new CoefficientFunction.
Note that all operations are done on the unit interval it is the 
reference time that is fixed. 
    """

    if not isinstance(time, Parameter):
      if isinstance(obj,GridFunction) or isinstance(obj,ProxyFunction):
        if time == 0:
          return obj.Operator("fix_tref_bottom")
        elif time == 1: 
          return obj.Operator("fix_tref_top")
        elif isinstance(obj,GridFunction):
          return fix_tref_gf(obj,time,*args,**kwargs)
      elif isinstance(obj,ProxyFunction):
        return fix_tref_proxy(obj,time,*args,**kwargs)

    if isinstance(obj,CoefficientFunction):
      return fix_tref_coef(obj,time,*args,**kwargs)
    else:
      raise Exception("obj is not a CoefficientFunction")

def fix_t_coef(obj,time,*args,**kwargs):
  print("WARNING: fix_t_coef is deprecated. Use \"fix_tref_coef\" instead. \n         Note that operators act w.r.t. the reference time intervals.")
  return fix_tref_coef(obj,time,*args,**kwargs)
def fix_t_gf(obj,time,*args,**kwargs):
  print("WARNING: fix_t_gf is deprecated. Use \"fix_tref_gf\" instead. \n         Note that operators act w.r.t. the reference time intervals.")
  return fix_tref_gf(obj,time,*args,**kwargs)
def fix_t_proxy(obj,time,*args,**kwargs):
  print("WARNING: fix_t_proxy is deprecated. Use \"fix_tref_proxy\" instead. \n         Note that operators act w.r.t. the reference time intervals.")
  return fix_tref_proxy(obj,time,*args,**kwargs)


def fix_t(obj,time,*args,**kwargs):
  """
  Deprecated: use "fix_tref" instead
  """
  print("WARNING: fix_t is deprecated. Use \"fix_tref\" instead. \n         Note that operators act w.r.t. the reference time intervals.")
  return fix_tref(obj,time,*args,**kwargs)

import ngsolve
from ngsolve.internal import *


class DummyScene:
  def __init__(self):
    pass
  def Redraw(self,blocking=False):
    ngsolve.Redraw(blocking=blocking)
dummy_scene = DummyScene()


def DrawDiscontinuous_std(StdDraw,levelset, fneg, fpos, *args, **kwargs):
    def StdDrawWithDummyScene(cf,*args,**kwargs):
        StdDraw(cf,*args,**kwargs)
        return dummy_scene

    if "deformation" in kwargs and StdDraw.__module__ == "ngsolve.solve":
        args2 = list(args[:])
        args2[1] = "deformation_"+args[1]
        StdDraw(kwargs["deformation"],*tuple(args2),**kwargs)
        visoptions.deformation=1
    if not "sd" in kwargs:
        kwargs["sd"] = 5
        
    return StdDrawWithDummyScene(IfPos(levelset,fpos,fneg),*args,**kwargs)
    
def DrawDiscontinuous_webgui(WebGuiDraw,levelset, fneg, fpos, *args, **kwargs):
    fneg = CoefficientFunction(fneg)
    fpos = CoefficientFunction(fpos)
    if fneg.dim > 1 or fpos.dim > 1:
        print("webgui discontinuous vis only for scalar functions a.t.m., switching to IfPos variant")
    else:
        return WebGuiDraw(CoefficientFunction((levelset,fpos,fneg,0)),eval_function="value.x>0.0?value.y:value.z",*args,**kwargs)
    return DrawDiscontinuous_std(WebGuiDraw,levelset, fneg, fpos, *args, **kwargs)

from functools import partial
def MakeDiscontinuousDraw(Draw):
    """
Generates a Draw-like visualization function. If Draw is from the webgui, a special evaluator is used to draw a pixel-sharp discontinuity otherwise an IfPos-CoefficientFunction is used.     
    """
    if (Draw.__module__ == "ngsolve.webgui"):
        ret = partial(DrawDiscontinuous_webgui,Draw)
    else:
        ret = partial(DrawDiscontinuous_std,Draw)
    ret.__doc__ ="""
    Visualization method for functions that are non-smooth across 
    level set interfaces. Effectively calls """+Draw.__module__+""".Draw with
    a few manipulations. 

    Parameters
    ----------
    levelset : CoefficientFunction
        (scalar) CoefficientFunction that describes the (implicit) geometry 
    fneg : CoefficientFunction
        CoefficientFunction that is evaluated where the level set function is negative
    fpos : CoefficientFunction
        CoefficientFunction that is evaluated where the level set function is positive
        
    deformation : deformation (optional)
        (vectorial) CoefficientFunction that describes the deformation of the background mesh
    *remainder* : *
        all remainder arguments are passed to """ +Draw.__module__ +".Draw"
    return ret

class NoDeformation:
    """
Dummy deformation class. Does nothing to the mesh. Has two dummy members:
  * lset_p1 : ngsolve.GridFunction
    The piecewise linear level set function 
  * deform : ngsolve.GridFunction
    A zero GridFunction (for compatibility with netgen Draw(... deformation=)) 
    """

    def __init__(self, mesh=None, levelset=None):
        self.deform = GridFunction(H1(mesh, order=1, dim=mesh.dim), "dummy_deform")
        self.deform.vec.data[:] = 0.0
        if levelset != None:
            if mesh == None:
                raise Exception("need mesh")
            self.lset_p1 = GridFunction(H1(mesh, order=1))
            InterpolateToP1(levelset, self.lset_p1)

    def __enter__(self):
        return self.lset_p1

    def __exit__(self, type, value, tb):
        pass


try:
    __IPYTHON__
    from ipywidgets import interact, FloatSlider
    from ngsolve.webgui import Draw
    def TimeSlider_Draw(cf,mesh,*args,**kwargs):
        ts = Parameter(0)
        if not isinstance(cf,CoefficientFunction):
            cf = CoefficientFunction(cf)
        scene = Draw(fix_tref(cf,ts),mesh,*args,**kwargs); 
        def UpdateTime(time): 
            ts.Set(time); scene.Redraw()
        return interact(UpdateTime,time=FloatSlider(description="tref:", 
                                                    continuous_update=False,
                                                    min=0,max=1,step=.025))
    def TimeSlider_DrawDC(cf1,cf2,cf3,mesh,*args,**kwargs):
        """
Draw a (reference) time-dependent function that is discontinuous across an 
interface described by a level set function. Change reference time through
widget slider.

        Args:
            cf1 (CoefficientFunction): level set function
            cf2 (CoefficientFunction): function to draw where lset is negative
            cf3 (CoefficientFunction): function to draw where lset is positive
            mesh (Mesh): Mesh

        Returns:
            widget element that allows to vary the reference time. 
        """
        DrawDC = MakeDiscontinuousDraw(Draw)
        if not isinstance(cf1,CoefficientFunction):
            cf1=CoefficientFunction(cf1)
        if not isinstance(cf2,CoefficientFunction):
            cf2=CoefficientFunction(cf2)
        if not isinstance(cf3,CoefficientFunction):
            cf3=CoefficientFunction(cf3)
        ts = Parameter(0)
        scene = DrawDC(fix_tref(cf1,ts),fix_tref(cf2,ts),fix_tref(cf3,ts),mesh,*args,**kwargs); 
        def UpdateTime(time): 
            ts.Set(time); scene.Redraw()
        return interact(UpdateTime,time=FloatSlider(description="tref:", 
                                                    continuous_update=False,
                                                    min=0,max=1,step=.025))
    import ngsolve.webgui
    DrawDC = MakeDiscontinuousDraw(ngsolve.webgui.Draw)
        
except:
    def TimeSlider_Draw(cf,mesh,*args,**kwargs):
      print("TimeSlider_Draw only available in ipython mode")
    def TimeSlider_DrawDC(cf1,cf2,cf3,mesh,*args,**kwargs):
      print("TimeSlider_DrawDC only available in ipython mode")
    import ngsolve
    DrawDC = MakeDiscontinuousDraw(ngsolve.Draw)


_dCut_raw = CutDifferentialSymbol(VOL)
_dFacetPatch_raw = FacetPatchDifferentialSymbol(VOL)

def dFacetPatch(**kwargs):
    """
    Differential symbol for facet patch integrators.

    Parameters
    ----------
    definedon : Region
        Domain description on where the integrator is defined.
    deformation : ngsolve.GridFunction
        Mesh deformation that is applied during integration. Default: None.
    definedonelements : ngsolve.BitArray
        Allows integration only on a set of facets
        that are marked True. Default: None.
    time_order : int
        Order in time that is used in the space-time integration.
        Default: time_order=-1 means that no space-time rule will be
        applied. This is only relevant for space-time discretizations.
    tref : double
        Turn spatial integration into space-time integration with 
        fixed time tref.

    Returns
    -------
      FacetPatchDifferentialSymbol(VOL)
    """
    if "element_vb" in kwargs or "element_boundary" in kwargs \
       or "skeleton" in kwargs:
        raise Exception("facet patch integrators are fixed to facet patches")
    return _dFacetPatch_raw(**kwargs)


def dCut(levelset, domain_type, order=None, subdivlvl=None, time_order=-1,
         levelset_domain=None, **kwargs):
    """
    Differential symbol for cut integration.

    Parameters
    ----------
    levelset : ngsolve.GridFunction
        The level set fct. describing the geometry 
        (desirable: P1 approximation).
    domain_type : {POS, IF, NEG, mlset.DomainTypeArray}
        The domain type of interest.
    order : int
        Modify the order of the integration rule used.
    subdivlvl : int
        Number of additional subdivision used on cut elements to
        generate the cut quadrature rule. Note: subdivlvl >0 only
        makes sense if you don't provide a P1 level set function
        and no isoparametric mapping is used.
    definedon : Region
        Domain description on where the integrator is defined.
    vb : {VOL, BND, BBND}
        Integration on mesh volume or its (B)boundary. Default: VOL
        (if combined with skeleton=True VOL refers to interior facets
                                        BND refers to boundary facets)
    element_boundary : bool
        Integration on each element boundary. Default: False
    element_vb : {VOL, BND, BBND}
        Integration on each element or its (B)boundary. Default: VOL
        (is overwritten by element_boundary if element_boundary 
        is True)
    skeleton : bool
        Integration over element-interface. Default: False.
    deformation : ngsolve.GridFunction
        Mesh deformation that is applied. Default: None.
    definedonelements : ngsolve.BitArray
        Allows integration only on elements or facets (if skeleton=True)
        that are marked True. Default: None.
    time_order : int
        Order in time that is used in the space-time integration.
        Default: time_order=-1 means that no space-time rule will be
        applied. This is only relevant for space-time discretizations.
    tref : float
        turns a spatial integral resulting in spatial integration rules
        into a space-time quadrature rule with fixed reference time tref
    levelset_domain : dict
        description of integration domain through a dictionary 
        (deprecated).

    Returns
    -------
        CutDifferentialSymbol(VOL)
    """
    if levelset_domain is not None and type(levelset_domain) == dict:
        lsetdom = levelset_domain
    else:
        lsetdom = {"levelset": levelset, "domain_type": domain_type}
    if order is not None and "order" not in lsetdom.keys():
        lsetdom["order"] = order
    if subdivlvl is not None and "subdivlvl" not in lsetdom.keys():
        lsetdom["subdivlvl"] = subdivlvl
    if time_order > -1 and "time_order" not in lsetdom.keys():
        lsetdom["time_order"] = time_order
    if "tref" in kwargs:
        lsetdom["tref"] = kwargs["tref"]
        del kwargs["tref"]

    return _dCut_raw(lsetdom, **kwargs)


def dxtref(mesh, order=None, time_order=-1, **kwargs):
    """
    Differential symbol for the integration over all elements extruded by
    the reference interval [0,1] to space-time prisms.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        The spatial mesh.
        The domain type of interest.
    order : int
        Modify the order of the integration rule used.
    definedon : Region
        Domain description on where the integrator is defined.
    vb : {VOL, BND, BBND}
        Integration on domains volume or boundary. Default: VOL
        (if combined with skeleton VOL means interior facets,
                                   BND means boundary facets)
    element_boundary : bool
        Integration on each element boundary. Default: False
    element_vb : {VOL, BND, BBND}
        Integration on each element or its (B)boundary. Default: VOL
        (is overwritten by element_boundary if element_boundary 
        is True)
    skeleton : bool
        Integration over element-interface. Default: False.
    deformation : ngsolve.GridFunction
        Mesh deformation. Default: None.
    definedonelements : ngsolve.BitArray
        Allows integration only on elements or facets (if skeleton=True)
        that are marked True. Default: None.
    time_order : int
        Order in time that is used in the space-time integration.
        Default: time_order=-1 means that no space-time rule will be
        applied. This is only relevant for space-time discretizations.

    Return
    ------
        CutDifferentialSymbol(VOL)
    """
    tFE = ScalarTimeFE(1)
    STFES = tFE*H1(mesh)
    gflset = GridFunction(STFES)  
    gflset.vec[:] = 1
    #for i in range(gflset.space.ndof):
    #    gflset.vec[i] = i+1


    lsetdom = {"levelset": gflset, "domain_type": POS}
    if order is not None:
        if type(order) != int:
            raise Exception("dxtref: order is not an integer! use keyword arguments for vb=VOL/BND.")
        lsetdom["order"] = order
    if type(time_order) != int:
        raise Exception("dxtref: time_order is not an integer! use keyword arguments for vb=VOL/BND.")
    if time_order > -1:
        lsetdom["time_order"] = time_order

    return _dCut_raw(lsetdom, **kwargs)

def dmesh(mesh=None,*args,**kwargs):
    """
    Differential symbol for the integration over all elements in the mesh.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        The spatial mesh.
        The domain type of interest.
    definedon : Region
        Domain description on where the integrator is defined.
    element_boundary : bool
        Integration on each element boundary. Default: False
    element_vb : {VOL, BND, BBND}
        Integration on each element or its (B)boundary. Default: VOL
        (is overwritten by element_boundary if element_boundary 
        is True)
    skeleton : bool
        Integration over element-interface. Default: False.
    deformation : ngsolve.GridFunction
        Mesh deformation. Default: None.
    definedonelements : ngsolve.BitArray
        Allows integration only on elements or facets (if skeleton=True)
        that are marked True. Default: None.
    tref : float
        turns a spatial integral resulting in spatial integration rules
        into a space-time quadrature rule with fixed reference time tref

    Return
    ------
        CutDifferentialSymbol(VOL)
    """
    if "tref" in kwargs:
        if mesh == None:
            raise Exception("dx(..,tref..) needs mesh")
        gflset = GridFunction(H1(mesh))
        gflset.vec[:] = 1
        lsetdom = {"levelset": gflset, "domain_type": POS, "tref" : kwargs["tref"]}
        if "order" in kwargs:
            lsetdom["order"] = kwargs["order"]
            del kwargs["order"]
        del kwargs["tref"]
        return _dCut_raw(lsetdom, **kwargs)
    else:
        return dx(*args,**kwargs)

def RestrictedBilinearForm(space=None,name="blf",element_restriction=None,facet_restriction=None,trialspace=None,testspace=None,**kwargs):
    """
    Creates a restricted bilinear form, which is bilinear form with a reduced MatrixGraph
    compared to the usual BilinearForm. BitArray(s) define on which elements/facets entries will be
    created.

    Use cases:

      * ghost penalty type stabilization:
        Facet-stabilization that are introduced only act on a few facets in the mesh. By providing the
        information on the corresponding facets, these additional couplings will only be introduced
        where necessary.

      * fictitious domain methods (domain decomposition methods):
        When PDE problems are only solved on a part of a domain while a finite element space is used
        that is still defined on the whole domain, a BitArray can be used to mark the 'active' part of
        the mesh. 

    Parameters

    space: ngsolve.FESpace
      finite element space on which the bilinear form is defined. If trial space and test space are different 
      they can be specified using the trialspace and testspace arguments.

    name : string
      name of the bilinear form

    element_restriction : ngsolve.BitArray
      BitArray defining the 'active mesh' element-wise

    facet_restriction : ngsolve.BitArray
      BitArray defining the 'active facets'. This is only relevant if FESpace has DG-terms (dgjumps=True)
    
    trialspace : ngsolve.FESpace
      finite element space on which the bilinear form is defined
      (trial space).
    
    testspace : ngsolve.FESpace
      finite element space on which the bilinear form is defined
      (test space).

    kwargs : keyword arguments
      kwargs are pasre to flags and passed to bilinearform 
    """

    argument_list = [name]
    if element_restriction != None:
        argument_list.append(element_restriction)
    if facet_restriction != None:
        argument_list.append(facet_restriction)

    if trialspace != None and testspace != None:
        if trialspace.is_complex and testspace.is_complex:
            return RestrictedBilinearFormComplex(trialspace,testspace,*argument_list,**kwargs)
        elif (not trialspace.is_complex) and (not testspace.is_complex):
            return RestrictedBilinearFormDouble(trialspace,testspace,*argument_list,**kwargs)
        else:
            raise Exception("Trialspace is real and testspace is complex or vice versa. RestrictedBilinearForm not implemented for this case!")
    else: 
        if space == None:
            raise Exception("No space given for RestrictedBilinearForm!")
        if space.is_complex:
            return RestrictedBilinearFormComplex(space,*argument_list,**kwargs)
        else:
            return RestrictedBilinearFormDouble(space,*argument_list,**kwargs)


# some global scope manipulations (monkey patches etc..):

# monkey patches
# print("|---------------------------------------------|")
# print("| ngsxfem applied monkey patches for          |")
# print("| GridFunction.Set, SymbolicLFI, SymbolicBFI. |")
# print("|---------------------------------------------|")
GridFunction.Set = SpaceTimeSet
SymbolicLFI = SymbolicLFIWrapper
SymbolicBFI = SymbolicBFIWrapper

# global scope definitions:
tref = ReferenceTimeVariable()
