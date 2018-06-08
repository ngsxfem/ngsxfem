from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *
from ngsolve.utils import *

from xfem import *

class LevelSetMeshAdaptation:
    """
Class to compute a proper mesh deformation to improve a piecewise (multi-) linear level set
approximation to obtain a higher order accurate approximation.

Computes a continuous function that describes a shift between points that are on the (P1 )
approximated level set function and its higher order accurate approximation. The transformation is
only applied on elements where a level value inside a certain interval (lower,upper) exists.

The result is a deform function (D) which is computed pointwise as

1)phi_lin( Psi(x) ) = phi_h(x)

  with Psi(x) = x + d(x) qn(x) =: D(x)

for all x on 'cut' elements

with

  phi_h : self.lset_ho
    the higher order level set function

  phi_lin : self.lset_p1
    the P1 level set function

  Psi : Id + self.deform
    the resulting deformation

  qn : self.qn
    normal direction field

This class holds its own members for the higher order an lower order (P1) approximation of the level
set function and only depends on the input of a mesh and a CoefficientFunction (the levelset
function) (and options).
    """

    order_deform = 2
    order_qn = 2
    order_lset = 2

    def __init__(self, mesh, order = 2, lset_lower_bound = 0, lset_upper_bound = 0, threshold = -1, discontinuous_qn = False, eps_perturbation = 1e-14, heapsize=1000000):
        """
The computed deformation depends on different options:

  order : int
    order of deformation GridFunction (ideally) order + 1 is the accuracy of the geometry
    approximation after applying the deformation on the mesh

  lset_lower_bound: float
    smallest relevant level set value to define the 'cut' elements where the mapping should be
    applied
    
  lset_upper_bound: float
    largest relevant level set value to define the 'cut' elements where the mapping should be
    applied
    
  threshold: float
    maximum (pointwise) value for d(x)/h in the mapping
      Psi(x) = x + d(x) qn(x)
    of the mesh transformation. A small value might be necessary if the geometry is only coarsely
    approximated to avoid irregular meshes after a corresponding mesh deformation.

  discontinuous_qn: boolean
    As an approximation for the normal direction we use n_h = nabla phi_h (gradient of higher order
    level set approximation) depending on the discontinuous_qn flag this normal field will be
    projected onto a continuous finite element space or not.

  eps_perturbation: float
    epsilon perturbation that is used to interpolate to P1

  heapsize : int
    heapsize for local computations.
        """
        self.order_deform = order
        self.order_qn = order
        self.order_lset = order

        self.lset_lower_bound = lset_lower_bound
        self.lset_upper_bound = lset_upper_bound
        self.threshold = threshold

        self.eps_perturbation = eps_perturbation
        
        self.v_ho = H1(mesh, order=self.order_lset)
        self.lset_ho = GridFunction (self.v_ho, "lset_ho")

        if (discontinuous_qn):
            self.v_qn = L2(mesh, order=self.order_qn, dim=mesh.dim)
        else:
            self.v_qn = H1(mesh, order=self.order_qn, dim=mesh.dim)
        self.qn = GridFunction(self.v_qn, "qn")
    
        self.v_p1 = H1(mesh, order=1)
        self.lset_p1 = GridFunction (self.v_p1, "lset_p1")

        self.v_def = H1(mesh, order=self.order_deform, dim=mesh.dim)
        self.deform = GridFunction(self.v_def, "deform")
        self.heapsize = heapsize

    def CalcDeformation(self, levelset, ba =None, blending=None):
        """
Compute the mesh deformation, s.t. isolines on cut elements of lset_p1 (the piecewise linear
approximation) are mapped towards the corresponding isolines of a given function

Parameters:

levelset : CoefficientFunction
  The coefficient function that prescribes where the piecewise linear iso lines should be mapped to.

blending : None/string/CoefficientFunction
  Option to apply the mesh deformation more localized on cut elements. Setting blending function to
  0 (CoefficientFunction(0.0) or None or "none") corresponds to applying the mapping on all points
  on cut elements completely. Using a blending function as a CoefficientFunction allows for a
  transition between the full application of the mapping (value 0) and no application of the mapping
  (value 1). There are predefined blending functions that are selectable via a string:
   * "quadratic":
     blending function that is 0 at the zero level set (of lset_p1) and increases quadratically with
     lset_p1. It is scaled with h, so that value 1 is not reached within cut elements.
   * "quartic":
     blending function that is 0 at the zero level set (of lset_p1) and increases like a fourth
     order polynomial with lset_p1. It is scaled with h, so that value 1 is not reached within cut
     elements.
        """
        self.v_ho.Update()
        self.lset_ho.Update()
        self.v_p1.Update()
        self.lset_p1.Update()
        self.v_qn.Update()
        self.qn.Update()
        self.v_def.Update()
        self.deform.Update()
        
        self.lset_ho.Set(levelset)
        self.qn.Set(self.lset_ho.Deriv())
        InterpolateToP1(self.lset_ho,self.lset_p1,eps_perturbation=self.eps_perturbation)
        if blending == None or blending == "none":
            blending = CoefficientFunction(0.0)
        elif blending == "quadratic":
            scale=sqrt(self.lset_p1.space.mesh.dim) * specialcf.mesh_size
            blending = self.lset_p1*self.lset_p1/( scale * scale)
        elif blending == "quartic":
            scale=sqrt(self.lset_p1.space.mesh.dim) * specialcf.mesh_size
            blending = self.lset_p1*self.lset_p1*self.lset_p1*self.lset_p1/(scale*scale*scale*scale)
            
        ProjectShift(self.lset_ho,
                     self.lset_p1,
                     self.deform,
                     self.qn,
                     ba,
                     blending,
                     lower=self.lset_lower_bound,
                     upper=self.lset_upper_bound,
                     threshold=self.threshold,
                     heapsize=self.heapsize)
        return self.deform


    # def CalcDistances(self, levelset,  lset_stats):
    #     """
    #     Compute largest distance
    #     """
    #     CalcDistances(levelset,self.lset_p1,self.deform,lset_stats)

    def CalcMaxDistance(self, levelset, heapsize=None):
        """
Compute approximated distance between of the isoparametrically obtained geometry.

See documentation of xfem.CalcMaxDistance 
        """
        if (heapsize == None):
            return CalcMaxDistance(levelset,self.lset_p1,self.deform,heapsize=self.heapsize)
        else:
            return CalcMaxDistance(levelset,self.lset_p1,self.deform,heapsize=heapsize)
        
    def MarkForRefinement(self, levelset = None, refine_threshold = 0.1, absolute = False):
        """
Marks elements for refinement where the geometry approximation is larger than a prescrbed relative
(or absolute) (to the mesh size) error. This will lead to a refinement in zones with high curvature.

Parameters:

  levelset : ngsolve.CoefficientFunction/None
    accurate level set function that is used to compute the error. If None self.lset_ho will be
    used.

  refine_threshold : float
    value that decides if a geometry error is large enough to mark for refinement or not

  absolute : boolean
    decides if the refine_threshold is an absolute value or if it is weighted with the mesh size
        """
        lset_stats = StatisticContainer()
        if levelset==None:
            levelset = self.lset_ho
        CalcDistances(lset_ho=levelset,lset_p1=self.lset_p1,deform=self.deform,stats=lset_stats,refine_threshold=refine_threshold,absolute=absolute)
        
#     def CalcDeformationError(self, deform_stats):
#         """
# Compute error between ideal and realized deformation
#         """
#         CalcDeformationError(self.lset_ho,self.lset_p1,self.deform,deform_stats,self.lset_lower_bound, self.lset_upper_bound)
        
