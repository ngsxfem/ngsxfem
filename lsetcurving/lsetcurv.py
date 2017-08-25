from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *
from ngsolve.utils import *

from xfem import *

class LevelSetMeshAdaptation:
    """
    Deformation from level set 
    """

    order_deform = 2
    order_qn = 2
    order_lset = 2

    def __init__(self, mesh, order = 2, lset_lower_bound = 0, lset_upper_bound = 0, threshold = -1, discontinuous_qn = False, heapsize=1000000):
        """
        Deformation
        """
        self.order_deform = order
        self.order_qn = order
        self.order_lset = order

        self.lset_lower_bound = lset_lower_bound
        self.lset_upper_bound = lset_upper_bound
        self.threshold = threshold
        
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

    def CalcDeformation(self, levelset, blending=None):
    """
    Compute the mesh deformation, s.t. isolines on cut elements of lset_p1 
    (the piecewise linear approximation) are mapped towards the corresponding
    isolines of a given function 
    
    Parameters:
    
    levelset : CoefficientFunction
      The coefficient function that prescribes where the piecewise linear iso 
      lines should be mapped to.
    
    blending : None/string/CoefficientFunction
      Option to apply the mesh deformation more localized on cut elements.
      Setting blending function to 0 (CoefficientFunction(0.0) or None or "none")
      corresponds to applying the mapping on all points on cut elements completely
      Using a blending function as a CoefficientFunction allows for a transition 
      between the full application of the mapping (value 0) and no application of 
      the mapping (value 1).
      There are predefined blending functions that are selectable via a string:
       * "quadratic":
         blending function that is 0 at the zero level set (of lset_p1) and 
         increases quadratically with lset_p1. It is scaled with h, so that value 
         1 is not reached within cut elements. 
       * "quartic":
         blending function that is 0 at the zero level set (of lset_p1) and 
         increases like a fourth order polynomial with lset_p1. It is scaled with 
         h, so that value 1 is not reached within cut elements. 
    """
        self.v_ho.Update()
        self.lset_ho.Update()
        self.v_p1.Update()
        self.lset_p1.Update()
        self.v_qn.Update()
        self.qn.Update()
        self.v_def.Update()
        self.deform.Update()
        
        self.lset_ho.Set(levelset,heapsize=self.heapsize)
        self.qn.Set(self.lset_ho.Deriv(),heapsize=self.heapsize)
        InterpolateToP1(self.lset_ho,self.lset_p1)
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
                     blending,
                     lower=self.lset_lower_bound,
                     upper=self.lset_upper_bound,
                     threshold=self.threshold,
                     heapsize=self.heapsize)
        return self.deform


    def CalcDistances(self, levelset,  lset_stats):
        """
        Compute largest distance
        """
        CalcDistances(levelset,self.lset_p1,self.deform,lset_stats)

    def CalcMaxDistance(self, levelset, heapsize=None):
        """
        Compute largest distance
        """
        if (heapsize == None):
            return CalcMaxDistance(levelset,self.lset_p1,self.deform,heapsize=self.heapsize)
        else:
            return CalcMaxDistance(levelset,self.lset_p1,self.deform,heapsize=heapsize)
        
    def MarkForRefinement(self, levelset, refine_threshold, absolute = False):
        """
        Compute largest distance
        """
        lset_stats = StatisticContainer()
        CalcDistances(lset_ho=levelset,lset_p1=self.lset_p1,deform=self.deform,stats=lset_stats,refine_threshold=refine_threshold,absolute=absolute)
        
    def CalcDeformationError(self, deform_stats):
        """
        Compute error in deformation
        """
        CalcDeformationError(self.lset_ho,self.lset_p1,self.deform,deform_stats,self.lset_lower_bound, self.lset_upper_bound)
        
