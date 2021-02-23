from ngsolve import *
from netgen.geom2d import SplineGeometry
from xfem import *
from math import pi
from xfem.lset_spacetime import *
ngsglobals.msg_level = 1

rect = SplineGeometry()
rect.AddRectangle([-0.6, -1], [0.6, 1])
ngmesh = rect.GenerateMesh(maxh=0.5, quad_dominated=False)
mesh = Mesh(ngmesh)

k_s = 2
k_t = 2
time_order = 2 * k_t

# spatial FESpace for solution
fes1 = H1(mesh, order=k_s)
# time finite elements (nodal!)
tfe = ScalarTimeFE(k_t) 
tfe_i = ScalarTimeFE(k_t, skip_first_node=True) # interior
tfe_e = ScalarTimeFE(k_t, only_first_node=True) # exterior (inital values)
tfe_t = ScalarTimeFE(k_t-1)                     # test
# space-time finite element space
st_fes = SpaceTimeFESpace(fes1,tfe, flags = {"dgjumps": True})
st_fes_i = SpaceTimeFESpace(fes1,tfe_i, flags = {"dgjumps": True})
st_fes_e = SpaceTimeFESpace(fes1,tfe_e, flags = {"dgjumps": True})
st_fes_t = SpaceTimeFESpace(fes1,tfe_t, flags = {"dgjumps": True})

u_e = st_fes_e.TrialFunction()
v_t = st_fes_t.TestFunction()

lsetadap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space=2,
                                            order_time=2,
                                            threshold=0.5,
                                            discontinuous_qn=True)

delta_t =  1./8
told = Parameter(0.)
t = told + delta_t * tref

# level set geometry
# radius of disk (the geometry)
R = 0.5
# position shift of the geometry in time
rho = (1 / (pi)) * sin(2 * pi * t)
# convection velocity:
w = CoefficientFunction((0, rho.Diff(t)))
# level set
r = sqrt(x**2 + (y - rho)**2)
levelset = r - R

dfm = lsetadap_st.CalcDeformation(levelset,tref)
ci = CutInfo(mesh, time_order=time_order)
ci.Update(lsetadap_st.levelsetp1[INTERVAL], time_order=time_order)

dQ = delta_t * dCut(lsetadap_st.levelsetp1[INTERVAL], NEG, time_order=2 * k_t,
                    deformation=lsetadap_st.deformation[INTERVAL],
                    definedonelements=ci.GetElementsOfType(HASNEG))

a_e = BilinearForm(trialspace = st_fes_e, testspace = st_fes_t, check_unused=False, symmetric=False)
a_e += (InnerProduct(grad(u_e), grad(v_t))) * dQ

a_e.Assemble()
print(a_e.mat)

mesh.SetDeformation(lsetadap_st.deformation[INTERVAL])
a_e.Assemble()
print(a_e.mat)
