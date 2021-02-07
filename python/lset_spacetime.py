from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from numpy import pi


class LevelSetMeshAdaptation_Spacetime:
    """
    Deformation from level set in space-time based on the 
    LevelSetMeshAdaptation class.
    """

    order_deform = 2
    order_qn = 2
    order_lset = 2

    def __init__(self, mesh, order_space = 2, order_time = 1, lset_lower_bound = 0,
                 lset_upper_bound = 0, threshold = -1, discontinuous_qn = False, heapsize=1000000,periodic=False):
        """
        Deformation
        """
        self.order_deform = order_space
        self.order_qn = order_space
        self.order_lset = order_space
        self.order_time = order_time

        self.lset_lower_bound = lset_lower_bound
        self.lset_upper_bound = lset_upper_bound
        self.threshold = threshold
        self.periodic = periodic
        
        self.mesh = mesh
        self.v_ho = H1(mesh, order=self.order_lset)
        self.lset_ho_node = GridFunction (self.v_ho, "lset_ho_node")
        self.ndof_node = len(self.lset_ho_node.vec)

        if (discontinuous_qn):
            self.v_qn = L2(mesh, order=self.order_qn, dim=mesh.dim)
        else:
            self.v_qn = H1(mesh, order=self.order_qn, dim=mesh.dim)
        self.qn = GridFunction(self.v_qn, "qn")
    
        self.v_p1 = H1(mesh, order=1)
        self.lset_p1_node = GridFunction (self.v_p1, "lset_p1_node")
        self.ndof_node_p1 = len(self.lset_p1_node.vec)

        if self.periodic:
            self.v_def = Periodic(H1(mesh, order=self.order_deform, dim=mesh.dim))
        else:
            self.v_def = H1(mesh, order=self.order_deform, dim=mesh.dim)
        
        self.deform_node = GridFunction(self.v_def, "deform_node")
        self.heapsize = heapsize
        
        # Spacetime
        self.tfe = ScalarTimeFE(self.order_time) 
        
        self.v_ho_st = SpaceTimeFESpace(self.v_ho,self.tfe)
        self.lset_ho = GridFunction(self.v_ho_st)
        
        self.v_p1_st = SpaceTimeFESpace(self.v_p1,self.tfe)
        self.lset_p1 = GridFunction(self.v_p1_st)
        
        self.v_def_st = SpaceTimeFESpace(self.v_def,self.tfe)
        self.deform = GridFunction(self.v_def_st)
        
        self.ci = CutInfo(mesh)
        
        self.hasneg_spacetime = BitArray(self.ci.GetElementsOfType(NEG))
        self.hasneg_spacetime[:] = False
        self.haspos_spacetime = BitArray(self.ci.GetElementsOfType(POS))
        self.haspos_spacetime[:] = False
        self.hasif_spacetime = BitArray(self.ci.GetElementsOfType(IF))
        self.hasif_spacetime[:] = False
        
        self.v_kappa_node = L2(mesh,order=0)
        self.v_kappa = SpaceTimeFESpace(self.v_kappa_node,self.tfe)
        self.kappa = GridFunction(self.v_kappa, "kappa")        
        
    def interpol_ho(self,levelset):
        times = [xi for xi in self.v_ho_st.TimeFE_nodes()]
        for i,ti in enumerate(times):
            self.lset_ho_node.Set(fix_t(levelset,ti))
            self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node].data = self.lset_ho_node.vec[:]

    def interpol_p1(self):
        for i in range(self.order_time + 1):
            self.lset_ho_node.vec[:].data = self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node]
            InterpolateToP1(self.lset_ho_node,self.lset_p1_node)
            self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1].data = self.lset_p1_node.vec[:]
            
    def CalcDeformation(self, levelset,calc_kappa = False):
        """
        Compute the deformation
        """
        self.v_ho.Update()
        self.lset_ho_node.Update()
        self.v_p1.Update()
        self.lset_p1_node.Update()
        self.v_qn.Update()
        self.qn.Update()
        self.v_def.Update()
        self.deform_node.Update()
        self.v_kappa_node.Update()
        self.v_kappa.Update()
        self.kappa.Update()
        
        self.interpol_ho(levelset)
        self.interpol_p1()
                
        for i in  range(len(self.v_ho_st.TimeFE_nodes())):
            self.lset_p1_node.vec[:].data = self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1]
            self.ci.Update(self.lset_p1_node)
            if calc_kappa:
                self.kappa.vec[i*self.v_kappa_node.ndof : (i+1)*self.v_kappa_node.ndof].data = self.ci.GetCutRatios(VOL)
        
        
        self.ci.Update(self.lset_p1,time_order=self.order_time)
        self.hasneg_spacetime[:] = False
        self.hasneg_spacetime |= self.ci.GetElementsOfType(NEG)
        self.haspos_spacetime[:] = False
        self.haspos_spacetime |= self.ci.GetElementsOfType(POS)
        self.hasif_spacetime[:] = False
        self.hasif_spacetime |= self.ci.GetElementsOfType(IF)
        
        
        for i in range(self.order_time + 1):
            self.lset_ho_node.vec[:].data = self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node]
            self.qn.Set(self.lset_ho_node.Deriv())
            self.lset_p1_node.vec[:].data = self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1]
            ProjectShift(self.lset_ho_node, self.lset_p1_node, self.deform_node, self.qn, self.hasif_spacetime, None, self.lset_lower_bound, 
                         self.lset_upper_bound, self.threshold, heapsize=self.heapsize)
            self.deform.vec[i*self.ndof_node : (i+1)*self.ndof_node].data = self.deform_node.vec[:]
        return self.deform
            

    def CalcMaxDistance(self, levelset, order=None, time_order=None, heapsize=None):
        """
Compute approximated distance between of the isoparametrically obtained geometry
(should be called in deformed state)
        """
        if order == None:
          order = 2 * self.order_qn
        if time_order == None:
          time_order = 2 * self.order_time
        if heapsize == None:
          heapsize = self.heapsize
        lset_dom = {"levelset": self.lset_p1, "domain_type" : IF, "order": order, "time_order" : time_order}
        minv, maxv = IntegrationPointExtrema(lset_dom, self.mesh, levelset, heapsize=heapsize)
        return max(abs(minv),abs(maxv))