from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

from ngsolve.internal import *
from xfem import *
from math import pi


class LevelSetMeshAdaptation_Spacetime:
    """
Class to compute a proper mesh deformation to improve a piecewise (multi-) 
linear (in space!) level set approximation to obtain a higher order 
accurate approximation.

Computes a continuous function that describes a shift between points that 
are on the (P1-in-space ) approximated level set function and its higher 
order accurate approximation. The transformation is only applied on elements 
where a level value inside a certain interval (lower,upper) exists.

The result is a space-time finite element deformation (`deform`). For each 
fixed time the behavior is as for an `LevelSetMeshAdaptation` object, i.e.


1)phi_lin( Psi(x) ) = phi_h(x)

  with Psi(x) = x + d(x) qn(x) =: D(x)

for all x on 'cut' elements

with

  phi_h : self.lset_ho
    the higher order space-time level set function

  phi_lin : self.lset_p1
    the P1-in-space space-time level set function

  Psi : Id + self.deform
    the resulting deformation

  qn : self.qn
    normal direction field

This class holds its own members for the higher order and lower order
 (P1-in-space) space-time approximation of the level set function and 
 only depends on the input of a mesh and a CoefficientFunction 
 (the levelset function) (and options).

 A LevelSetMeshAdaptation_Spacetime can also be used as a context
 manager. In this case, the mesh deformation is applied on the mesh
 inside the context.
    """

    order_deform = 2
    order_qn = 2
    order_lset = 2


    @TimeFunction
    def __init__(self, mesh, order_space = 2, order_time = 1, lset_lower_bound = 0,
                 lset_upper_bound = 0, threshold = -1, discontinuous_qn = False, heapsize=1000000,periodic=False):

        """
The computed deformation depends on different options:

  order_space : int
    order of deformation GridFunction (ideally) order + 1 is the accuracy of the geometry
    approximation after applying the deformation on the mesh

  order_time : int
    order in time for space-time tensor-product spaces

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

  levelset : CoefficientFunction or None(default)
    If a level set function is prescribed the deformation is computed right away. Otherwise the 
    computation is triggered only at calls for `CalcDeformation`.
        """

        self.gf_to_project = []
        self.gf_to_project_tmp = [] # list for temporary copies
        
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
        self.lset_p1_top = CreateTimeRestrictedGF(self.lset_p1,1.0)
        self.lset_p1_bottom = CreateTimeRestrictedGF(self.lset_p1,0.0)        
        
        self.v_def_st = SpaceTimeFESpace(self.v_def,self.tfe)
        self.deform = GridFunction(self.v_def_st)
        self.deform_top = CreateTimeRestrictedGF(self.deform,1.0)
        self.deform_bottom = CreateTimeRestrictedGF(self.deform,0.0)
        
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

        self.deform_last_top = CreateTimeRestrictedGF(self.deform,1.0)

        class MeshDeformationContext:
            def __init__(self,lsetadap_st,type="std"):
                self.lsetadap_st = lsetadap_st
                self.type = type
            def __enter__(self):
                if self.type == "std":
                    self.lsetadap_st.mesh.SetDeformation(self.lsetadap_st.deform)
                    return self.lsetadap_st.lset_p1
                elif self.type == "top":   
                    self.lsetadap_st.mesh.SetDeformation(self.lsetadap_st.deform_top)
                    return self.lsetadap_st.lset_p1_top
                elif self.type == "bottom":   
                    self.lsetadap_st.mesh.SetDeformation(self.lsetadap_st.deform_bottom)
                    return self.lsetadap_st.lset_p1_bottom
                else:
                    raise Exception("no suitable deform type given")
            def __exit__(self, type, value, tb):
                self.lsetadap_st.mesh.UnsetDeformation()
        self.top = MeshDeformationContext(self,"top")
        self.bottom = MeshDeformationContext(self,"bottom")

        self.levelsetp1 = {INTERVAL : self.lset_p1, BOTTOM : self.lset_p1_bottom, TOP : self.lset_p1_top}
        self.deformation = {INTERVAL : self.deform, BOTTOM : self.deform_bottom, TOP : self.deform_top}

        
    def ProjectOnUpdate(self,gf,update_domain=None):
        """
When the LevelsetMeshAdaptation class generates a new deformation (due to 
a new level set function) all GridFunction that have been stored through
`ProjectOnUpdate` will be projected from the previous to the new mesh
by an essentially local projection (Oswald projection of the shifted 
evaluation)

  gfu : GridFunction (or list of GridFunctions)
    GridFunction(s) to store for later deformation updates.
        """      
        if isinstance(gf,list):
            self.gf_to_project.extend([(gfa, update_domain) for gfa in gf])
            self.gf_to_project_tmp.extend([GridFunction(gfa.space) for gfa in gf])
        else:
            self.gf_to_project.append((gf, update_domain))
            self.gf_to_project_tmp.append(GridFunction(gf.space))

    @TimeFunction
    def ProjectGFs(self):
        """
ProjectGFs projects all stored GridFunctions to the currect deformation.
This function is typically only called in the `CalcDeformation` unless
`CalcDeformation` is called with the argument `dont_project_gfs = True`.
        """      

        for (gf, update_domain),gftmp in zip(self.gf_to_project, self.gf_to_project_tmp):
            gftmp.vec.data = gf.vec
            if update_domain is None:
                if (gf.space.mesh.ne != self.deform_bottom.space.mesh.ne):
                    raise Exception("different meshes...")
                gf.Set(shifted_eval(gftmp, back = self.deform_last_top, forth = self.deform_bottom))
            else:
                gf.Set(shifted_eval(gftmp, back = self.deform_last_top, forth = self.deform_bottom), definedonelements=update_domain)
            #print("updated ", gf.name)

    @TimeFunction
    def interpol_ho(self,levelset):
        """
Internal function that is called inside `CalcDeformation`. It projects 
the space-time coefficient function to a space-time tensor product finite 
element function.
        """
        times = [xi for xi in self.v_ho_st.TimeFE_nodes()]
        for i,ti in enumerate(times):
            self.lset_ho_node.Set(fix_tref(levelset,ti))
            self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node].data = self.lset_ho_node.vec[:]

    @TimeFunction
    def interpol_p1(self):
        """
Internal function that is called inside `CalcDeformation`. It projects 
the space-time tensor product finite element approximation lset_ho 
to a piecewise linear-in-space approximation through interpolation
(in space).
        """
        for i in range(self.order_time + 1):
            self.lset_ho_node.vec[:].data = self.lset_ho.vec[i*self.ndof_node : (i+1)*self.ndof_node]
            InterpolateToP1(self.lset_ho_node,self.lset_p1_node)
            self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1].data = self.lset_p1_node.vec[:]

    @TimeFunction
    def CalcDeformation(self, levelset, calc_kappa = False, dont_project_gfs = False):
        """
Compute the space-time mesh deformation, s.t. isolines on cut elements of lset_p1 (the piecewise 
linear-in-space tensor product approximation) are mapped towards the corresponding isolines of 
a given space-time function

Parameters:

levelset : CoefficientFunction
  The coefficient function that prescribes where the piecewise linear-in-space iso lines should 
  be mapped to.

calc_kappa : bool
  Compute the cut ratio of a spatial element as a space-time function. Populates `self.kappa`.

dont_project_gfs : bool
  Usually, all stored GridFunctions are projected on the new deformed mesh. This flag
  can be used to avoid the call of the projection step (`ProjectGFs`).
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
        self.deform_last_top.Update()               
        self.interpol_ho(levelset)
        self.interpol_p1()
        RestrictGFInTime(spacetime_gf=self.deform,reference_time=1.0,space_gf=self.deform_last_top)   


        for i in  range(len(self.v_ho_st.TimeFE_nodes())):
            self.lset_p1_node.vec[:].data = self.lset_p1.vec[i*self.ndof_node_p1 : (i+1)*self.ndof_node_p1]
            if calc_kappa:
                self.ci.Update(self.lset_p1_node)
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


        RestrictGFInTime(spacetime_gf=self.lset_p1,reference_time=0.0,space_gf=self.lset_p1_bottom)
        RestrictGFInTime(spacetime_gf=self.lset_p1,reference_time=1.0,space_gf=self.lset_p1_top)
        RestrictGFInTime(spacetime_gf=self.deform,reference_time=0.0,space_gf=self.deform_bottom)   
        RestrictGFInTime(spacetime_gf=self.deform,reference_time=1.0,space_gf=self.deform_top)   

        if not dont_project_gfs:
            self.ProjectGFs()
        return self.deform
            
    def __enter__(self):
      self.mesh.SetDeformation(self.deform)
      return self.lset_p1

    def __exit__(self, type, value, tb):
      self.mesh.UnsetDeformation()

    def levelset_domain(self, domain_type = IF, time_type = INTERVAL):
        """
Return a levelset_domain dictionary.
      Args:
          domain_type (DOMAIN_TYPE: NEG/POS/IF, optional): domain type for level set doamin. Defaults to IF.
          time_type (BOTTOM/TOP/INTERVAL): type of integral on the space-time domain

      Returns:
          dict: levelset domain
        """        
        return { "levelset" : self.levelsetp1[time_type], "domain_type" : domain_type}



    def Integrator(self, SymbolicFI, domain_type, time_type, form, time_order = None, definedonelements = None):
        """
Convenience function to construct space-time Symbolic(Cut)LFI/Symbolic(Cut)BFI

        Args:
            SymbolicFI (SymbolicLFI/SymbolicBFI): type of integrator (LFI vs. BFI)
            domain_type (DOMAIN_TYPE: NEG/POS/IF): domain type to integrate on
            time_type (BOTTOM/TOP/INTERVAL): type of integral on the space-time domain
            form (CoefficientFunction): integrand
            time_order (int, optional): integration order in time. Defaults to None.
            defineonelement (BitArray, optional): elements to define the linear form on. Defaults to None.
        """    
        if time_order == None:
            time_order = 2 * self.order_time
        
        fi = SymbolicFI(levelset_domain = self.levelset_domain(domain_type,time_type), 
                         form = form, 
                         time_order=time_order, 
                         deformation=self.deformation[time_type])
        if definedonelements != None:
            fi.SetDefinedOnElements(definedonelements)       
        return fi

    def LFI(self, domain_type, time_type, form, time_order = None, definedonelements = None):
        """
Convenience function to construct Symbolic(Cut)LFI based on the domain_type 
and the known space-time level set function.

        Args:
            domain_type (DOMAIN_TYPE: NEG/POS/IF): domain type to integrate on
            time_type (BOTTOM/TOP/INTERVAL): type of integral on the space-time domain
            form (CoefficientFunction): integrand
            time_order (int, optional): integration order in time. Defaults to None.
            defineonelement (BitArray, optional): elements to define the linear form on. Defaults to None.
        """    
        return self.Integrator(SymbolicLFI, domain_type, time_type, form, time_order, definedonelements)

    def BFI(self, domain_type, time_type, form, time_order = None, definedonelements = None):
        """
Convenience function to construct Symbolic(Cut)LFI based on the domain_type 
and the known space-time level set function.

        Args:
            domain_type (DOMAIN_TYPE: NEG/POS/IF): domain type to integrate on
            time_type (BOTTOM/TOP/INTERVAL): type of integral on the space-time domain
            form (CoefficientFunction): integrand
            time_order (int, optional): integration order in time. Defaults to None.
            defineonelement (BitArray, optional): elements to define the bilinear form on. Defaults to None.
        """    
        return self.Integrator(SymbolicBFI, domain_type, time_type, form, time_order, definedonelements)

    @TimeFunction
    def Integrate(self, domain_type, time_type, cf, order = 5, time_order = None):
        """
Convenience function to Integrate on cut space-time domain.

        Args:
            domain_type (DOMAIN_TYPE: NEG/POS/IF): domain type to integrate on
            time_type (BOTTOM/TOP/INTERVAL): type of integral on the space-time domain
            cf (CoefficientFunction): integrand
            order (int, optional): integration order in space. Defaults to 5.
            time_order (int, optional): integration order in time. Defaults to None.

        Returns:
            float: return value of (cut) integral
        """
        if time_order == None:
            time_order = 2 * self.order_time
        self.mesh.SetDeformation(self.deformation[time_type])
        fi = Integrate(levelset_domain = self.levelset_domain(domain_type,time_type), 
                       mesh = self.mesh, cf = cf, order = order, 
                       time_order=time_order)
        self.mesh.UnsetDeformation()
        return fi


    @TimeFunction
    def CalcMaxDistance(self, levelset, order=None, time_order=None, heapsize=None):
        """
Compute approximated distance between of the isoparametrically obtained geometry
(should be called in deformed state)

        Args:
            levelset (CoefficientFunction): implicit geometry description, 
              assumed to be a signed distance function.
            order (int, optional): integration order in space. Defaults to 5.
            time_order (int, optional): integration order in time. Defaults to None.
            heapsize (int, optional): heap size for local memory. Defaults to None.
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
