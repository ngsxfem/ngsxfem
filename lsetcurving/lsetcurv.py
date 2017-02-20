from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *

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

    def CalcDeformation(self, levelset):
        """
        Compute the deformation
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
        ProjectShift(self.lset_ho, self.lset_p1, self.deform, self.qn, self.lset_lower_bound, self.lset_upper_bound, self.threshold, heapsize=self.heapsize);
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
        
