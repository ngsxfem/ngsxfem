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
        
    def interpol_ho(self,levelset,t): #,tstart,delta_t):
        #times = [tstart + delta_t * xi for xi in self.v_ho_st.TimeFE_nodes()]
        times = [xi for xi in self.v_ho_st.TimeFE_nodes()]
        for i,ti in enumerate(times):
            #t.Set(ti)
            t.FixTime(ti)
            #print("i, ti: ", i, ti)
            self.lset_ho_node.Set(levelset)
            self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node].data = self.lset_ho_node.vec[:]
        #t.Set(tstart)
        t.UnfixTime()

    def interpol_p1(self):
        for i in range(self.order_time + 1):
            self.lset_ho_node.vec[:].data = self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node]
            InterpolateToP1(self.lset_ho_node,self.lset_p1_node)
            self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1].data = self.lset_p1_node.vec[:]
            
    def CalcDeformation(self, levelset,t,calc_kappa = False):
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
        
        self.interpol_ho(levelset,t) #,tstart,delta_t)
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
            
    def CalcMaxDistance(self, levelset,t, given_pts = []):
        """
        Compute largest distance
        """
        if given_pts:
            time_quad = given_pts
        else:
            time_quad = self.v_ho_st.TimeFE_nodes() 
        # times = [tstart + delta_t * xi for xi in time_quad]
        max_dists = []
        # for ti,xi in zip(times,time_quad):
        for xi in time_quad:
        
            RestrictGFInTime(self.lset_p1,xi,self.lset_p1_node)
            RestrictGFInTime(self.deform,xi,self.deform_node)
            t.FixTime(xi)
            # t.Set(ti) 
            self.v_def_st.SetTime(xi)
            self.v_ho_st.SetTime(xi)
            max_dists.append(CalcMaxDistance(levelset,self.lset_p1_node,self.deform_node,heapsize=self.heapsize))
            #max_dists.append(CalcMaxDistance(self.lset_ho,self.lset_p1,self.deform,heapsize=self.heapsize))
        t.UnfixTime()
        self.v_def_st.SetOverrideTime(False)
        self.v_ho_st.SetOverrideTime(False)
        # t.Set(tstart)
        return max(max_dists)
      
       
## geometry        
#square = SplineGeometry()
#square.AddRectangle([0,0],[2,2],bc=1)
#ngmesh = square.GenerateMesh(maxh=0.03, quad_dominated=False)
#mesh = Mesh (ngmesh)
#
## data
#t = Parameter(0)
#lset = CoefficientFunction( sqrt( (x-1-0.25*sin(pi*t))*(x-1-0.25*sin(pi*t))+(y-1)*(y-1)) - 0.5  )
#
#
#
#def SolveProblem(mesh,delta_t,k_s=2,k_t=1):
#    max_nodes = []
#    max_interm = []
#    tstart = 0
#    tnew = 0
#    tend = 1
#    t.Set(tnew)
#    lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = k_t,
#                                             threshold=0.5, discontinuous_qn=True)
#    
#    while tend - tnew > delta_t/2:
#        t.Set(tnew)
#        dfm = lset_adap_st.CalcDeformation(lset,t,tnew,delta_t,calc_kappa=True) 
#        Draw(lset_adap_st.kappa,mesh,"kappa")
#        max_nodes.append(lset_adap_st.CalcMaxDistance(lset,t,tnew,delta_t))
#        max_interm.append(lset_adap_st.CalcMaxDistance(
#                              lset,t,tnew,delta_t,[i*0.1 for i in range(11)]))
#        tnew += delta_t 
#    return max(max_nodes),max(max_interm)
#
#
#def StudyConvergence(delta_t=0.5,max_rfs=2,where = "space"):   
#    max_dist_nodes = []
#    max_dist_interm = []
#    ref_lvl = 0
#    while ref_lvl <= max_rfs:  
#        
#        if where == "space" and ref_lvl > 0:
#            mesh.Refine()
#        elif where == "time" and ref_lvl > 0:
#            delta_t = delta_t / 2  
#        e1,e2 = SolveProblem(mesh,delta_t,k_s=3,k_t=2)
#        max_dist_nodes.append(e1)
#        max_dist_interm.append(e2)
#        ref_lvl = ref_lvl + 1
#    print("Studying convergence in: " + where)
#    print("Max-Dist nodes = {0}".format(max_dist_nodes))
#    if min(max_dist_nodes):
#        eoc_nodes = [ log(max_dist_nodes[i-1]/max_dist_nodes[i])/log(2) for i in range(1,len(max_dist_nodes))]
#        print("Eoc nodes = {0}".format(eoc_nodes))
#    print("Max-Dist intermediate = {0}".format(max_dist_interm))
#    if min(max_dist_interm):
#        eoc_interm = [ log(max_dist_interm[i-1]/max_dist_interm[i])/log(2) for i in range(1,len(max_dist_interm))]
#        print("Eoc intermediate = {0}".format(eoc_interm))
#            
#            
#StudyConvergence(delta_t=0.01,max_rfs=3,where = "time")
#        
#                                            
#lset_adap_st.interpol_ho(lset,t,tstart,delta_t)
#lset_ho = lset_adap_st.lset_ho
#lset_adap_st.interpol_p1()
#lset_p1 = lset_adap_st.lset_p1
#dfm = lset_adap_st.CalcDeformation(lset,t,tstart,delta_t)
#print("Max-Dist nodes = {0}".format(lset_adap_st.CalcMaxDistance(lset,t,tstart,delta_t)))
#print("Max-Dist intermediate points = {0}".format(lset_adap_st.CalcMaxDistance(
#                              lset,t,tstart,delta_t,[i*0.1 for i in range(11)])))
#
## Plotting
#visoptions.deformation = 1
#lset_adap_st.v_ho_st.SetTime(0.0)    
#lset_adap_st.v_p1_st.SetTime(0.0) 
#lset_adap_st.v_def_st.SetTime(0.0) 
#Draw(lset_ho,mesh,"lsetHO")
#Draw(lset_p1,mesh,"lsetP1")
#Draw(dfm,mesh,"deformation")
#input("")
#lset_adap_st.v_ho_st.SetTime(1.0) 
#lset_adap_st.v_p1_st.SetTime(1.0)  
#lset_adap_st.v_def_st.SetTime(1.0)  
#Redraw()   
