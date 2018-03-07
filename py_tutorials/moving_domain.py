"""
In this example we solve a moving domain problem with homogeneous Neumann boundary conditions.
The domain's evolution is described by a given level set function. We work with a simple background 
mesh that is unfitted to the moving domain. The problem is treated by a fictitious domain approach
combined with a Space-Time-DG method. To obtain higher order accuracy we employ a mesh trans-
formation at fixed points in time.

    domain: 
    -------  
    The background domain is chosen as [-1,-0.75]x[1,1.5] while the moving domain is described
    as the set of points where the level set function takes negative values.

    PDE problem:
    ------------  
    domain equations:
    u_t - u_xx - u_yy + w_1 * u_x + w_2 * u_y = f in moving domain
    boundary conditions:
                                        du/dn = 0 on boundary of moving domain

    discretization:
    ---------------
    Variational formulation:
    We use a Space-Time-DG method for the variational formulation. In order to reduce 
    computational complexity the problem is divided into time intervals (time slabs). Using
    a suited weak formulation on these time slabs allows to solve the problem time slab by
    time slab. This yields the computational structure of a time stepping scheme.

    Finite element space:
    The functions in our finite element space are linear combinations of tensor products
    between purely time respectively space dependent functions. The polynomial degrees in 
    space and time can be prescribed independently. 
    
    implementation:
    ---------------
    The integrals that appear in the variational formulation are first approximated
    by using quadrature in time, cf. [1]. This approach requires additional stabilization in 
    order to extend the discrete solution. To this end we apply Ghost-Penalty stabilization.
    
    linear systems:
    ---------------
    A (sparse) direct solver is applied to solve the arising linear systems.

    further information:
    ---------------
    A more detailed explanation of the implementation is given in a
    juypter-notebook available from [2].
    
    literature:
    -----------
    [1]: P. Hansbo, M Larson, and S. Zahedi. A CutFEM for coupled bulk-surface problems 
         on time-dependent domains. CMAME, 2016
    [2]: jupyter – Tutorials for ngsxfem. http://www.github.com/ngsxfem/ngsxfem-jupyter

"""

from ngsolve import *
from ngsolve.comp import *
from ngsolve.utils import *
from xfem import *
from netgen.geom2d import SplineGeometry
from xfem.lsetcurv import *
from math import pi


class quad_rule:

    def __init__(self,name,npoints):
        '''Constructor of quadrature rule. 
           This class is used for approximating the
           time integrals.
        
        Parameters
        ----------
        name : str
            Name of the quadrature rule.
        npoints : int
            Number of quadrature points.
        '''
        self.name = name
        self.npoints = npoints

        # available quadrature rules
        gauss_radau = {
            3: ([-1, (1-sqrt(6))/5, (1+sqrt(6))/5],
                [2/9, (16+sqrt(6))/18,(16-sqrt(6))/18]),
            4: ([-1, -0.575319, 0.181066, 0.822824],
                [0.125,0.657689,0.776387,0.440924])}
                
       # you may add your favourite quadrature rule here 
                
        if name == "Gauss-Radau":
            self.points = gauss_radau[npoints][0]
            self.weights = gauss_radau[npoints][1]

    def shifted_pts(self,a,b):
        '''Transformation of quadrature points to interval [a,b].
        
        Parameters
        ----------
        a : float
            Left end of the interval.
        b : float
            Right end of the interval.

        Returns
        -------
        list
            Transformed quadrature points.       
        '''
        if self.name == "Gauss-Radau":
            return [0.5*(b-a) * pt + 0.5*(b+a)  for pt in self.points]

    def shifted_weights(self,a,b):
        '''Transformation of quadrature weights to interval [a,b].
        
        Parameters
        ----------
        a : float
            Left end of the interval.
        b : float
            Right end of the interval.

        Returns
        -------
        list
            Transformed quadrature weights.
        '''
        if self.name == "Gauss-Radau":
            return [0.5*(b-a)*w for w in self.weights]

def dnjump(u,order,comp):
    '''Jump of normal derivative over element interfaces.
    
    Parameters
    ----------
    u : ProxyFunction
        Test or Trialfunction.
    order : int
            Order of the normal derivative
    comp : int
           Component of ProxyFunction.
           
    Returns
    -------
    ProxyFunction
        The jump of the normal derivative of the required order.    
    '''
    if order%2==0:
        return dn(u,order,comp) - dn(u.Other(),order,comp)
    else:
        return dn(u,order,comp) + dn(u.Other(),order,comp)

def power(u,p):
    if p == 0:
        return 1
    else:
        return u * power(u,p-1)

# geometry
square = SplineGeometry()
square.AddRectangle([-1,-0.75],[1,1.5],bc=1)
maxh = 0.1
mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=False))

# define quadrature rule
qr= quad_rule("Gauss-Radau",4)

# constants
r0 = 0.5
r1 = 0.5*pi/r0
tend  = 0.5

# the time
t = Parameter(0)

# for evaluation of time basis functions defined on reference interval [0,1]
tref = Parameter(0) 

# available basis functions in time and their derivatives
# you may expand this dictionary with your own basis 
basis_in_time = { "modal": ([lambda t: 1.0, lambda t: t, lambda t: (2*t-1)*(2*t-1)],[lambda t: 0.0, lambda t: 1, lambda t: 4*(2*t-1)]) }

# choice of basis                  
phi,d_phi_ref = basis_in_time["modal"]

rho =  CoefficientFunction((1/(pi))*sin(2*pi*t))
d_rho = CoefficientFunction(2*cos(2*pi*t))
w = CoefficientFunction((0,d_rho)) # convection

# the levelset function that describes the domain's evolution
levelset= CoefficientFunction(sqrt(x*x+(y-rho)*(y-rho)) -r0)

# right hand side
coeff_f = CoefficientFunction(  -(pi/r0)*r1*( sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*sin(r1*sqrt(x*x+(y-rho)*(y-rho))) - cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*cos(r1*sqrt(x*x+(y-rho)*(y-rho))) )  + (pi/r0)*cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*sin(r1*sqrt(x*x+(y-rho)*(y-rho)))*(1/sqrt(x*x+(y-rho)*(y-rho))) )

# exact solution for monitoring the error
exact = CoefficientFunction(cos(r1*sqrt(x*x+(y-rho)*(y-rho)))*cos(r1*sqrt(x*x+(y-rho)*(y-rho))) )

# initial condition
coef_u0 = cos(r1*sqrt(x*x+y*y))*cos(r1*sqrt(x*x+y*y))

# order of FE-Spaces
k_s = 2 #space
k_t = 2 # time

# FE-Spaces
V=H1(mesh, order = k_s, dirichlet=[],dgjumps = True)
u0 = GridFunction(V)
W = FESpace([V for l in range(k_t+1)],flags = {"dgjumps" : True}) # added DG Jumps for Ghost-penalty terms
gfu = GridFunction(W)
Draw(u0,mesh,"unew")

# stabilization parameter for Ghost-penalty
gamma_stab = [10**-l for l in range(1,k_s+1)]



class SlabDiscretization:
    '''
    Class for discretization on the time slab.
    '''
    def Update(self, delta_t):
        '''Sets up the problem on the time slabs.
    
        Parameters
        ----------
        delta_t : float
            The time step.
        '''
        self.active_domain_key = "hasneg_spacetime"
        self.bfs = []
        self.lfs = []
        # for recording on which part of the mesh the integrators are defined
        self.marked_integrators = []

        self.points_time_ref = []
        self.weights_time_ref = []

        # add dummy time point (with weight 0) if initial time is not part of time quad:
        if not 0.0 in qr.shifted_pts(0,1):
            self.points_time_ref.append(0)
            self.weights_time_ref.append(0)
            
        self.points_time_ref.extend(qr.shifted_pts(0,1))
        self.weights_time_ref.extend(qr.shifted_weights(0,1))
    
        if not 1.0 in qr.shifted_pts(0,1):
           self.points_time_ref.append(1.0)
           self.weights_time_ref.append(0.0)
        print("points_time_ref", self.points_time_ref)
        print("weights_time_ref", self.weights_time_ref)

        # mesh adaptation
        self.lsetmeshadap = LevelSetMeshAdaptation(mesh, order=k_s, threshold=0.1, discontinuous_qn=True)
        self.levelset_domain = { "levelset" : self.lsetmeshadap.lset_p1, "domain_type" : NEG }
    
        # define Bi/Linearforms:
        u = W.TrialFunction()
        v = W.TestFunction()

        h = specialcf.mesh_size
        
        # define Trial/Testfunctions 
        # these are linear combinations of tensor products (time x space)
        u_ti =   sum( [ u[n]*phi[n](tref) for n in range(k_t +1 ) ]  )  
        v_ti =   sum( [ v[n]*phi[n](tref) for n in range(k_t +1 ) ]  )
        # time-derivative (integration by parts in time)
        d_v_ti = sum( [ v[n]*d_phi_ref[n](tref) for n in range(k_t +1 ) ]  ) 
        gradu_ti = CoefficientFunction((0,0))
        gradv_ti = CoefficientFunction((0,0))
        for n in range(k_t+1):
            gradu_ti += grad(u[n])*phi[n](tref)
            gradv_ti += grad(v[n])*phi[n](tref)
        # jump of the normal derivative (needed for Ghost-Penalty stabilization)
        u_dnjump_ti = [ sum( [ phi[n](tref)*dnjump(u[n],l+1,n) for n in range(k_t+1)  ]) for l in range(k_s) ]
        v_dnjump_ti = [ sum( [ phi[n](tref)*dnjump(v[n],l+1,n) for n in range(k_t+1)  ]) for l in range(k_s) ]
        
        # loop over discrete time points to set up the integrators       
        # the actual linear system is assembled in the SolveProblem function below. 
        for t_i, omega_i in zip(self.points_time_ref, self.weights_time_ref):
                           
            ai = BilinearForm(W,symmetric=False,check_unused=False)
            fi = LinearForm(W)
            
            marked_integrators_at_ti = {}
            marked_integrators_at_ti["hasneg_at_ti"] = []
            marked_integrators_at_ti["gpfacets_spacetime"] = []
            
            # inner terms
            inner_term =  CoefficientFunction(0)
            inner_term_rhs = CoefficientFunction(0)
            if omega_i != 0.0:
                # diffusion term
                inner_term += delta_t*omega_i*gradu_ti*gradv_ti
                # convection
                inner_term += -delta_t*omega_i*u_ti*w*gradv_ti
                # time-derivative
                inner_term += -omega_i*u_ti*d_v_ti
                # right hand side
                inner_term_rhs += delta_t*omega_i*coeff_f*v_ti
            
            # boundary terms at the time slab's bottom/top 
            # these stem from integration by parts of the time derivative            
            if t_i == 0.0:
                inner_term_rhs += u0 * v_ti                   
            if t_i == 1.0:
                inner_term += u_ti*v_ti
                
            form_inner = SymbolicBFI(self.levelset_domain,form = inner_term)
            form_inner_rhs = SymbolicLFI(self.levelset_domain, form=inner_term_rhs)
            
            ai += form_inner
            fi += form_inner_rhs
            marked_integrators_at_ti["hasneg_at_ti"].append(form_inner)
            marked_integrators_at_ti["hasneg_at_ti"].append(form_inner_rhs)
                            
            if omega_i != 0.0:
                # Ghost-Penalty Terms
                for l in range(k_s):
                    gp_term = delta_t * omega_i * gamma_stab[l] * power(h,2*l+1) * u_dnjump_ti[l] * v_dnjump_ti[l]
                    form_new_gp = SymbolicBFI( gp_term, skeleton=True )
                    ai += form_new_gp
                    marked_integrators_at_ti["gpfacets_spacetime"].append(form_new_gp)
                    
            self.bfs.append(ai)
            self.lfs.append(fi)
            self.marked_integrators.append(marked_integrators_at_ti)
            
    def __init__(self):
        pass
            
        
            

def SolveProblem(sd, delta_t ):
    '''This function carries out the time stepping.
    
        Parameters
        ----------
        sd : SlabDiscretization
             Provides the integrators on the time slabs.
        delta_t : float
            Size of the time step.
            
        Returns
        -------
        float
            The L2-Error at the final point in time. 
    '''

    sd.Update(delta_t)
    u0.Update()
    gfu.Update() # also updates W

    # keeping track of time
    tstart  = 0
    told = tstart
    tnew = told
    t.Set(tstart)
    
    deformation = sd.lsetmeshadap.CalcDeformation(levelset)
    mesh.SetDeformation(deformation)            
    u0.Set(coef_u0) # set initial condition
    mesh.UnsetDeformation()
    
    # Dummy-Assemble
    a =  BilinearForm(W,symmetric=False)
    f = LinearForm(W)
    a.Assemble()
    f.Assemble()
    amat = a.mat.CreateMatrix()
    fvec = f.vec.CreateVector()

    # provides information about how the elements are cut (unfitted method)
    ci = CutInfo(mesh)

    # loop for the time stepping scheme
    while tend - tnew > delta_t/2:

        #updates
        tnew +=  delta_t

        # clear storage
        amat.AsVector()[:] = 0
        fvec[:] = 0

        # marker types:
        markers = {}
        # BitArrays for marking an element's type w.r.t. the partition
        # of the domain provided by the piecewise linear levelset function
        hasneg_spacetime = BitArray(ci.GetElementsOfType(NEG))
        hasneg_spacetime[:] = False
        haspos_spacetime = BitArray(ci.GetElementsOfType(POS))
        haspos_spacetime[:] = False
        hasif_spacetime = BitArray(ci.GetElementsOfType(IF))
        hasif_spacetime[:]= False
        
        # Loop over time points ti.
        # Elements are classified as hasneg_spacetime if they have 
        # been part of the active mesh on at least one of the time 
        # points ti. The linear system for the time slab will involve 
        # the Dofs belonging to these elements.
        for t_i in sd.points_time_ref:
            t.Set(told + t_i * delta_t)
            deformation = sd.lsetmeshadap.CalcDeformation(levelset)
            ci.Update(sd.lsetmeshadap.lset_p1)
            hasneg_spacetime |= ci.GetElementsOfType(HASNEG)
            haspos_spacetime |= ci.GetElementsOfType(HASPOS)
            hasif_spacetime |= ci.GetElementsOfType(IF)

        # elements that have been in haspos and in hasneg should be marked as hasif!
        jumpels = BitArray(hasneg_spacetime)
        jumpels &= haspos_spacetime
        hasif_spacetime |= jumpels

        markers["hasneg_spacetime"] = hasneg_spacetime
        markers["haspos_spacetime"] = haspos_spacetime
        markers["hasif_spacetime"] = hasif_spacetime

        # Draw(BitArrayCF(hasif_spacetime),mesh,"hasifOrjump")
        # Draw(BitArrayCF(hasneg_spacetime),mesh,"hasneg")

        # determine the facets on which Ghost-Penalty stabilization is applied
        ba_facets = GetFacetsWithNeighborTypes(mesh,a=hasneg_spacetime,b=hasif_spacetime,bnd_val_a=False,bnd_val_b=False)
        markers["gpfacets_spacetime"] = ba_facets
        
        # We are now on a fixed time slab and start approximating the integrals 
        # by first applying quadrature in time.
        for i, t_i in enumerate(sd.points_time_ref):

            tref.Set(t_i)
            t.Set(told + t_i * delta_t)
            
            # calculate the mesh deformation at the current point in time
            deformation = sd.lsetmeshadap.CalcDeformation(levelset)
            
            # collect information about the current cut-situation
            ci.Update(sd.lsetmeshadap.lset_p1)
            hasneg_ti = BitArray(ci.GetElementsOfType(HASNEG))
            haspos_ti = BitArray(ci.GetElementsOfType(HASPOS))

            markers["hasneg_at_ti"] = hasneg_ti
            markers["haspos_at_ti"] = haspos_ti
            markers["hasif_at_ti"] = ci.GetElementsOfType(IF) 


            # Collect the integrators from slab discretization and inform them
            # about the current cut-situation.
            for key in markers:
                if key in sd.marked_integrators[i] and len(sd.marked_integrators[i][key])>0:
                    for form in sd.marked_integrators[i][key]:
                        form.SetDefinedOnElements(markers[key])

            # Assemble the matrices at time point ti and add them 
            # to the total system.
            # mesh adaptation
            mesh.SetDeformation(deformation)            
            sd.bfs[i].Assemble()
            sd.lfs[i].Assemble()            
            # unset mesh adaptation
            mesh.UnsetDeformation()
            amat.AsVector().data += sd.bfs[i].mat.AsVector()
            fvec.data += sd.lfs[i].vec

        # solve linear system
        active_dofs = GetDofsOfElements(W,markers[sd.active_domain_key])
        gfu.vec.data = amat.Inverse(active_dofs,"umfpack") * fvec
        tref.Set(1)
        # solution at the time slab's top
        u0.vec[:] = 0.0
        for l in range(k_t+1):
            u0.vec.data += phi[l](1.0)*gfu.components[l].vec
            
        Redraw(blocking=True)

        # update
        told = tnew
 
    # measure the error
    t.Set(tend)
    # mesh adaptation
    deformation = sd.lsetmeshadap.CalcDeformation(levelset)
    mesh.SetDeformation(deformation)
    l2error = sqrt(Integrate(sd.levelset_domain,(u0-exact)*(u0-exact),mesh))    
    print("Time: {0}, Error: {1}".format(tnew,l2error))
    # unset mesh adaptation
    mesh.UnsetDeformation()
    return l2error

# construct slab discretization
sd = SlabDiscretization()

# solve the problem with a time step of delta_t
with TaskManager():
    SolveProblem(sd, delta_t = 0.01)
