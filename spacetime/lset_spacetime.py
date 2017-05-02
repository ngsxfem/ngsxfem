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
                 lset_upper_bound = 0, threshold = -1, discontinuous_qn = False, heapsize=1000000):
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
        
    def interpol_ho(self,levelset,t,tstart,delta_t):
        times = [tstart + delta_t * xi for xi in self.v_ho_st.TimeFE_nodes().NumPy()]
        for i,ti in enumerate(times):
            t.Set(ti)
            self.lset_ho_node.Set(levelset)
            self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node] = self.lset_ho_node.vec[:] 

    
    def interpol_p1(self):
        for i in range(self.order_time + 1):
            self.lset_ho_node.vec[:] = self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node]
            self.lset_p1_node.Set(self.lset_ho_node)
            self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1] = self.lset_p1_node.vec[:]
            
    def CalcDeformation(self, levelset,t,tstart,delta_t):
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
        
        self.interpol_ho(levelset,t,tstart,delta_t)
        self.interpol_p1()
        for i in range(self.order_time + 1):
            self.lset_ho_node.vec[:] = self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node]
            self.qn.Set(self.lset_ho_node.Deriv(),heapsize=self.heapsize)
            self.lset_p1_node.vec[:] = self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1]
            ProjectShift(self.lset_ho_node, self.lset_p1_node, self.deform_node, self.qn, self.lset_lower_bound, 
                         self.lset_upper_bound, self.threshold, heapsize=self.heapsize)
            self.deform.vec[i*self.ndof_node : (i+1)*self.ndof_node] = self.deform_node.vec[:]
        return self.deform
            
    def CalcMaxDistance(self, levelset,t,tstart,delta_t):
        """
        Compute largest distance
        """
        times = [tstart + delta_t * xi for xi in self.v_ho_st.TimeFE_nodes().NumPy()]
        max_dists = []
        for i,ti in enumerate(times):
            t.Set(ti)
            self.lset_p1_node.vec[:] = self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1]           
            self.deform_node.vec[:] = self.deform.vec[i*self.ndof_node : (i+1)*self.ndof_node]
            max_dists.append(CalcMaxDistance(levelset,self.lset_p1_node,self.deform_node,heapsize=self.heapsize))
        return max(max_dists)
      
       
# geometry        
square = SplineGeometry()
square.AddRectangle([0,0],[2,2],bc=1)
ngmesh = square.GenerateMesh(maxh=0.1, quad_dominated=False)
mesh = Mesh (ngmesh)

# data
t = Parameter(0)
tstart = 0
delta_t = 0.5
lset = CoefficientFunction( sqrt( (x-1-0.25*t)*(x-1-0.25*t)+(y-1)*(y-1)) - 0.5  )

# Spacetime interpolation
ref_lvl = 0
for i in range(ref_lvl):
    mesh.Refine()
lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = 2, order_time = 1,
                                              threshold=0.1, discontinuous_qn=True)                                            

                                              
lset_adap_st.interpol_ho(lset,t,tstart,delta_t)
lset_ho = lset_adap_st.lset_ho
lset_adap_st.interpol_p1()
lset_p1 = lset_adap_st.lset_p1
dfm = lset_adap_st.CalcDeformation(lset,t,tstart,delta_t)
print("Max-Dist = {0}".format(lset_adap_st.CalcMaxDistance(lset,t,tstart,delta_t)))

# Plotting
visoptions.deformation = 1
lset_adap_st.v_ho_st.SetTime(0.0)    
lset_adap_st.v_p1_st.SetTime(0.0) 
lset_adap_st.v_def_st.SetTime(0.0) 
Draw(lset_ho,mesh,"lsetHO")
Draw(lset_p1,mesh,"lsetP1")
Draw(dfm,mesh,"deformation")
input("")
lset_adap_st.v_ho_st.SetTime(1.0) 
lset_adap_st.v_p1_st.SetTime(1.0)  
lset_adap_st.v_def_st.SetTime(1.0)  
Redraw()   
