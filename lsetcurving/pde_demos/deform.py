from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *

from libngsxfem_py.xfem import *

class LevelSetMeshAdaptation:
    """
    Deformation from level set 
    """

    order_deform = 2
    order_qn = 2
    order_lset = 2

    # lset
    # lset_ho
    # lset_p1

    # deform
    
    def __init__(self, mesh, order = 2, lset_lower_bound = 0, lset_upper_bound = 0, threshold = -1):
        """
        Deformation
        """
        self.order_deform = order
        self.order_qn = 1
        self.order_lset = order

        self.lset_lower_bound = lset_lower_bound
        self.lset_upper_bound = lset_upper_bound
        self.threshold = threshold
        
        self.v_ho = FESpace("h1ho", mesh, order=self.order_lset)
        self.lset_ho = GridFunction (self.v_ho, "lset_ho")

        self.v_qn = FESpace("h1ho", mesh, order=self.order_qn, flags = { "vec" : True })
        self.qn = GridFunction(self.v_qn, "qn")
    
        self.v_p1 = FESpace("h1ho", mesh, order=1)
        self.lset_p1 = GridFunction (self.v_p1, "lset_p1")

        self.v_def = FESpace("h1ho", mesh, order=self.order_deform, flags = { "vec" : True })
        self.deform = GridFunction(self.v_def, "deform")
        

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
        
        self.lset_ho.Set(levelset)
        self.qn.Set(self.lset_ho.Deriv())
        InterpolateToP1(self.lset_ho,self.lset_p1)
        ProjectShift(self.lset_ho, self.lset_p1, self.deform, self.qn, self.lset_lower_bound, self.lset_upper_bound, self.threshold);
        return self.deform


    def CalcDistances(self, levelset,  lset_stats):
        """
        Compute largest distance
        """
        CalcDistances(levelset,self.lset_p1,self.deform,lset_stats)

    def CalcMaxDistance(self, levelset):
        """
        Compute largest distance
        """
        return CalcMaxDistance(levelset,self.lset_p1,self.deform)
        
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
        
