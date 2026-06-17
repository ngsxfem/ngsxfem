from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *
from ngsolve.utils import *

from xfem import *

# `from xfem import *` shadows ngsolve's Integrate with the cut-integration one;
# keep a handle to the plain ngsolve Integrate for the boundary defect measure.
from ngsolve.comp import Integrate as _ngs_integrate

class LevelSetMeshAdaptation:
    """
Class to compute a proper mesh deformation to improve a piecewise (multi-) linear level set
approximation to obtain a higher order accurate approximation.

Computes a continuous function that describes a shift between points that are on the (P1 )
approximated level set function and its higher order accurate approximation. The transformation is
only applied on elements where a level value inside a certain interval (lower,upper) exists.

The result is a deform function (D) which is computed pointwise as

1)phi_lin( x ) = phi_h(Psi(x))

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

A LevelSetMeshAdaptation can also be used as a context  manager. In this case, 
the mesh deformation is applied on the mesh inside the context.
    """

    order_deform = 2
    order_qn = 2
    order_lset = 2

    @TimeFunction
    def __init__(self, mesh, order = 2, lset_lower_bound = 0.0, lset_upper_bound = 0.0,
                 threshold = -1.0, discontinuous_qn = False,
                 eps_perturbation = 1e-14, heapsize=1000000,
                 levelset = None,
                 boundary_tangential = False, boundary_normal = None):
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

  levelset : CoefficientFunction or None(default)
    If a level set function is prescribed the deformation is computed right away.

  boundary_tangential : bool / str / list of str / Region (default False)
    Handle the case where the zero level set *crosses the domain boundary*: the
    search direction (qn) and the deformation (u) are made tangential to the
    selected boundaries (u.n = 0), so the deformed boundary stays on the
    (possibly curved) domain boundary and the higher-order cut accuracy is kept.
    Selects on which boundaries: False (off); True / ".*" (all); a name / regex
    ("left|bottom"); a list (["left", "bottom"]); or a boundary ``Region``
    (``mesh.Boundaries(...)``).  Implies an L2 qn field.

    Assumes the boundary is smooth where the interface crosses, and needs a mesh
    with a sound boundary parametrisation for ``specialcf.normal`` (a netgen
    mesh; ``MakeStructured2DMesh`` has degenerate left/right edges).

  boundary_normal : CoefficientFunction or None (default)
    Unused (kept for backward compatibility); the normal is ``specialcf.normal``.
        """
        self.order_deform = order
        self.order_qn = order
        self.order_lset = order

        self.lset_lower_bound = lset_lower_bound
        self.lset_upper_bound = lset_upper_bound
        self.threshold = threshold

        self.eps_perturbation = eps_perturbation

        self.boundary_tangential = boundary_tangential
        self.boundary_normal = boundary_normal
        # the qn tangential correction is applied via an element-local L2 Set
        # (whole mesh), which leaves the interior unchanged only for an L2 qn
        if boundary_tangential and not discontinuous_qn:
            discontinuous_qn = True

        self.mesh = mesh
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
        self.deform_last = GridFunction(self.v_def, "deform_last")
        self.heapsize = heapsize
        self.gf_to_project = []

        # boundary-tangential bookkeeping (built lazily, rebuilt on dof change)
        self._bt_normals = None       # [(Region, nrm, nc, full-facet els)]
        self._bt_def_region = None    # Region of selected boundaries
        self._bt_qn_tmp = None        # reusable scratch on qn / deform spaces
        self._bt_def_tmp = None
        self._bt_ndof = None
        if boundary_tangential:
            self._BuildBoundaryProjections()

        if levelset != None:
          self.CalcDeformation(levelset)

    # ----------------------------------------------------------------------- #
    #  boundary-tangential handling (for level sets crossing the boundary)
    #
    #  Per selected boundary we store, via specialcf.normal:
    #    * nrm : a VectorH1 field with the normal on the boundary, used on the
    #            boundary trace for the deformation (u.n = 0);
    #    * nc  : its piecewise-constant (L2 order 0) per-element projection, and
    #    * els : the volume elements with a full facet on that boundary;
    #  the search direction qn is made tangential on els with the *constant* nc
    #  (an exact, per-element projection).  Using the varying/normalised nrm in
    #  the volume, or correcting vertex-touching elements, over-corrects qn and
    #  folds the mesh on coarse grids -- hence nc on full-facet elements only.
    # ----------------------------------------------------------------------- #
    def _BoundaryPattern(self):
        bt = self.boundary_tangential
        if bt is True or bt == ".*":
            return ".*"
        return bt if isinstance(bt, str) else "|".join(bt)

    def _BuildBoundaryProjections(self):
        """Per selected boundary: normal field (nrm), its per-element constant
        projection (nc) and the full-facet volume elements (els)."""
        mesh = self.mesh
        D = mesh.dim
        bt = self.boundary_tangential
        self._bt_def_region = bt if isinstance(bt, Region) \
            else mesh.Boundaries(self._BoundaryPattern())
        bnames = mesh.GetBoundaries()
        sel = sorted({bnames[i] for i, on in enumerate(self._bt_def_region.Mask()) if on})
        els = {nm: BitArray(mesh.ne) for nm in sel}
        for b in els.values():
            b.Clear()
        for be in mesh.Elements(BND):
            nm = bnames[be.index]
            if nm in els:
                for e in mesh[be.elementnode].elements:
                    els[nm][e.nr] = True
        self._bt_normals = []
        for nm in sel:
            region = mesh.Boundaries(nm)
            nrm = GridFunction(VectorH1(mesh, order=self.order_deform))
            nrm.Set(specialcf.normal(D), definedon=region)
            nc = GridFunction(L2(mesh, order=0, dim=D))
            nc.Set(nrm)
            self._bt_normals.append((region, nrm, nc, els[nm]))
        # reusable scratch GridFunctions for the additive corrections
        self._bt_qn_tmp = GridFunction(self.qn.space)
        self._bt_def_tmp = GridFunction(self.deform.space)
        self._bt_ndof = self.v_def.ndof

    def _BoundaryProjectionsValid(self):
        return self._bt_ndof == self.v_def.ndof

    def _ApplyBoundaryTangentialQn(self):
        # qn <- qn - (qn.n) n on the full-facet boundary elements, n = constant
        # per-element normal nc (exact projection: (q.nc) nc / (nc.nc))
        if not self._BoundaryProjectionsValid():
            self._BuildBoundaryProjections()
        corr = self._bt_qn_tmp
        for _, nrm, nc, els in self._bt_normals:
            corr.Set(-(self.qn * nc) * nc / (nc * nc + 1e-30), definedonelements=els)
            self.qn.vec.data += corr.vec

    def _ApplyBoundaryTangentialDeform(self):
        # u <- u - (u.n) n on each boundary trace (nrm is unit there) -> u.n = 0
        if not self._BoundaryProjectionsValid():
            self._BuildBoundaryProjections()
        corr = self._bt_def_tmp
        for region, nrm, nc, els in self._bt_normals:
            corr.Set(-(self.deform * nrm) * nrm, definedon=region)
            self.deform.vec.data += corr.vec

    def CalcBoundaryNormalDefect(self):
        """L2 norm of u.n over the treated boundaries (~0 if u stays on them)."""
        if not self.boundary_tangential:
            return 0.0
        if not self._BoundaryProjectionsValid():
            self._BuildBoundaryProjections()
        n = specialcf.normal(self.mesh.dim)
        return sqrt(_ngs_integrate((self.deform * n) ** 2, self.mesh,
                                   definedon=self._bt_def_region))

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
        if isinstance(gf,list) or isinstance(gf,tuple):
            self.gf_to_project.extend([(gfa, update_domain) for gfa in gf])
        else:
            self.gf_to_project.append((gf, update_domain))

    @TimeFunction
    def ProjectGFs(self):
        """
ProjectGFs projects all stored GridFunctions to the currect deformation.
This function is typically only called in the `CalcDeformation` unless
`CalcDeformation` is called with the argument `dont_project_gfs = True`.
        """      

        timer = Timer("ProjectGFs")
        timer.Start()
        for (gf, update_domain) in self.gf_to_project:
            # make tmp copy 
            gfcopy = GridFunction(gf.space)
            gfcopy.vec.data = gf.vec
            if update_domain is None:
                gf.Set(shifted_eval(gfcopy, back = self.deform_last, forth = self.deform))
            else:
                gf.Set(shifted_eval(gfcopy, back = self.deform_last, forth = self.deform), definedonelements=update_domain)
            #print("updated ", gf.name)
        timer.Stop()

    @TimeFunction
    def CalcDeformation(self, levelset, ba =None, blending=None, dont_project_gfs=False):
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
dont_project_gfs : bool
  Usually, all stored GridFunctions are projected on the new deformed mesh. This flag
  can be used to avoid the call of the projection step (`ProjectGFs`).
        """
        self.v_ho.Update()
        self.lset_ho.Update()
        self.v_p1.Update()
        self.lset_p1.Update()
        self.v_qn.Update()
        self.qn.Update()
        self.v_def.Update()
        self.deform.Update()
        self.deform_last.Update()
        
        self.lset_ho.Set(levelset)
        self.qn.Set(self.lset_ho.Deriv())
        # (boundary) project the search direction onto the boundary tangent so
        # the deformation stays tangential where the interface crosses the
        # boundary instead of pushing the boundary off the domain
        if self.boundary_tangential:
            self._ApplyBoundaryTangentialQn()
        InterpolateToP1(self.lset_ho,self.lset_p1,eps_perturbation=self.eps_perturbation)
        if blending == None or blending == "none":
            blending = CoefficientFunction(0.0)
        elif blending == "quadratic":
            scale=sqrt(self.lset_p1.space.mesh.dim) * specialcf.mesh_size
            blending = self.lset_p1*self.lset_p1/( scale * scale)
        elif blending == "quartic":
            scale=sqrt(self.lset_p1.space.mesh.dim) * specialcf.mesh_size
            blending = self.lset_p1*self.lset_p1*self.lset_p1*self.lset_p1/(scale*scale*scale*scale)

        self.deform_last.vec.data = self.deform.vec
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

        # (boundary) remove the boundary-normal component of the deformation on
        # the boundary facets (u . n = 0): exact for straight facets, higher
        # order for curved ones
        if self.boundary_tangential:
            self._ApplyBoundaryTangentialDeform()

        if not dont_project_gfs:
            self.ProjectGFs()

        return self.deform

    def __enter__(self):
      self.mesh.SetDeformation(self.deform)
      return self.lset_p1

    def __exit__(self, type, value, tb):
      self.mesh.UnsetDeformation()

    @TimeFunction
    def CalcMaxDistance(self, levelset, order=None, heapsize=None, deform=False):
        """
Compute approximated distance between of the isoparametrically obtained geometry
(should be called in deformed state)
        """
        if order == None:
          order = 2 * self.order_qn
        if heapsize == None:
          heapsize = self.heapsize
        lset_dom = {"levelset": self.lset_p1, "domain_type" : IF, "order": order}
        if self.mesh.deformation or not deform:
          minv, maxv = IntegrationPointExtrema(lset_dom, self.mesh, levelset, heapsize=heapsize)
        else:
          self.mesh.deformation = self.deform
          minv, maxv = IntegrationPointExtrema(lset_dom, self.mesh, levelset, heapsize=heapsize)
          self.mesh.deformation = None
        return max(abs(minv),abs(maxv))

    def levelset_domain(self, domain_type = IF):
        """
  Return a levelset_domain dictionary.
        Args:
            domain_type (DOMAIN_TYPE, optional): domain type for level set doamin. Defaults to IF.

        Returns:
            dict: levelset domain
        """
        return { "levelset" : self.lset_p1, "domain_type" : domain_type}

    def Integrator(self, SymbolicFI, domain_type, form, definedonelements = None):
        """
  Convenience function to construct Symbolic(Cut)LFI/Symbolic(Cut)BFI
        """    
        fi = SymbolicFI(levelset_domain = self.levelset_domain(domain_type), 
                         form = form, 
                         deformation=self.deform)
        if definedonelements != None:
            fi.SetDefinedOnElements(definedonelements)       
        return fi

    def LFI(self, domain_type, form, definedonelements = None):
        """
  Convenience function to construct Symbolic(Cut)LFI based on the domain_type 
  and the known level set function.
        """    
        return self.Integrator(SymbolicLFI, domain_type, form, definedonelements)

    def BFI(self, domain_type, form, definedonelements = None):
        """
  Convenience function to construct Symbolic(Cut)BFI based on the domain_type
  and the known level set function.
        """    
        return self.Integrator(SymbolicBFI, domain_type, form, definedonelements)

    @TimeFunction
    def Integrate(self, domain_type, cf, order = 5):
        """
  Convenience function to Integrate on cut domain.
        """    
        self.mesh.SetDeformation(self.deform)
        fi = Integrate(levelset_domain = self.levelset_domain(domain_type), 
                       mesh = self.mesh, cf = cf, order = order)
        self.mesh.UnsetDeformation()
        return fi


    @TimeFunction
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
        
