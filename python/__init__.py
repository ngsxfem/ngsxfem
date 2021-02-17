"""
(ngs)xfem
=========

A module for unfitted discretizations in NGSolve

Modules:
xfem.lsetcurving ... isoparametric unfitted FEM
xfem.mlset ... multiple level sets
"""

from ngsolve import (L2, VOL, BitArray, CoefficientFunction, FESpace,
                     GridFunction, IfPos, LinearForm, Parameter)
from ngsolve.comp import Integrate as ngsolve_Integrate
from ngsolve.comp import ProxyFunction
from ngsolve.comp import SymbolicBFI as ngsolve_SymbolicBFI
from ngsolve.comp import SymbolicLFI as ngsolve_SymbolicLFI
from xfem.ngsxfem_py import *


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

def dt(func):
    """
Evaluates the time derivative of a Space-Time function.
    """
    add = func.Operator("dt")
    if add:
        return add
    raise Exception("cannot form dt")



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
    gfu_e_repl = GridFunction(space_fes)
    gfu_e_repl.Set( cf )
    gfu_e.vec[:].data = gfu_e_repl.vec

ngsolveSet = GridFunction.Set
def SpaceTimeSet(self, cf, *args, **kwargs):
    if (isinstance(self.space,CSpaceTimeFESpace)):
      cf = CoefficientFunction(cf)
      gfs = GridFunction(self.space.spaceFES)
      ndof_node = len(gfs.vec)
      j = 0
      for i,ti in enumerate(self.space.TimeFE_nodes()):
        if self.space.IsTimeNodeActive(i):
          ngsolveSet(gfs,fix_t(cf,ti), *args, **kwargs)
          self.vec[j*ndof_node : (j+1)*ndof_node].data = gfs.vec[:]
          j += 1
    else:
      ngsolveSet(self,cf, *args, **kwargs)

def fix_t(obj,time,*args,**kwargs):
    if not isinstance(time, Parameter):
      if isinstance(obj,GridFunction) or isinstance(obj,ProxyFunction):
        if time == 0:
          return obj.Operator("fix_t_bottom")
        elif time == 1: 
          return obj.Operator("fix_t_top")
        elif isinstance(obj,GridFunction):
          return fix_t_gf(obj,time,*args,**kwargs)
        elif isinstance(obj,ProxyFunction):
          return fix_t_proxy(obj,time,*args,**kwargs)

    if isinstance(obj,CoefficientFunction):
      return fix_t_coef(obj,time,*args,**kwargs)
    else:
      raise Exception("obj is not a CoefficientFunction")

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
    lsetp1 = None
    def __init__(self,mesh = None, levelset=None):
        if levelset != None:
            if mesh == None:
                raise Exception("need mesh");
            self.lsetp1 = GridFunction(H1(mesh))
            InterpolateToP1(levelset,self.lsetp1)

        pass
    def __enter__(self):
        return self.lsetp1
    def __exit__(self, type, value, tb):
        pass


try:
    __IPYTHON__
    from ipywidgets import FloatSlider, interact
    from ngsolve.webgui import Draw
    def TimeSlider_Draw(cf,mesh,*args,**kwargs):
        ts = Parameter(0)
        if not isinstance(cf,CoefficientFunction):
            cf = CoefficientFunction(cf)
        scene = Draw(fix_t(cf,ts),mesh,*args,**kwargs); 
        def UpdateTime(time): 
            ts.Set(time); scene.Redraw()
        return interact(UpdateTime,time=FloatSlider(description="tref:", 
                                                    continuous_update=True,
                                                    min=0,max=1,step=.025))
    def TimeSlider_DrawDC(cf1,cf2,cf3,mesh,*args,**kwargs):
        DrawDC = MakeDiscontinuousDraw(Draw)
        if not isinstance(cf1,CoefficientFunction):
            cf1=CoefficientFunction(cf1)
        if not isinstance(cf2,CoefficientFunction):
            cf2=CoefficientFunction(cf2)
        if not isinstance(cf3,CoefficientFunction):
            cf3=CoefficientFunction(cf3)
        ts = Parameter(0)
        scene = DrawDC(fix_t(cf1,ts),fix_t(cf2,ts),fix_t(cf3,ts),mesh,*args,**kwargs); 
        def UpdateTime(time): 
            ts.Set(time); scene.Redraw()
        return interact(UpdateTime,time=FloatSlider(description="tref:", 
                                                    continuous_update=True,
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

# some global scope manipulations (monkey patches etc..):

# monkey patches
GridFunction.Set = SpaceTimeSet
SymbolicLFI = SymbolicLFIWrapper
SymbolicBFI = SymbolicBFIWrapper

# global scope definitions:
tref = ReferenceTimeVariable()
