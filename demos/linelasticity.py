# %% ------------------------- LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

ngsglobals.msg_level = 2

# %% --------------------------- PARAMETERS ---------------------------------
# Domain corners
ll, ur = (0, -0.5), (2, 0.5)
# Initial mesh diameter
initial_maxh = 0.4
# Number of mesh bisections
nref = 4
# Order of finite element space
k = 1

# Stabilization parameter for ghost-penalty
gamma_s = 0.5
# Stabilization parameter for Nitsche
gamma_n = 10

# ----------------------------------- MAIN ------------------------------------
# %% Set up the level sets, exact solution and right-hand side


def level_sets():
    return [x - 1.5, y - 0.15, -0.15 - y]


nr_ls = len(level_sets())
rhs = CF((0, -0.1))

# %% Geometry and mesh
geo = SplineGeometry()
geo.AddRectangle(ll, ur, bcs=("bottom", "right", "top", "left"))
ngmesh = geo.GenerateMesh(maxh=initial_maxh)
for i in range(nref):
    ngmesh.Refine()
mesh = Mesh(ngmesh)


# %% Level set and cut-information
P1 = H1(mesh, order=1)
lsetsp1 = tuple(GridFunction(P1) for i in range(nr_ls))
for i, lsetp1 in enumerate(lsetsp1):
    InterpolateToP1(level_sets()[i], lsetp1)
    Draw(lsetp1, mesh, "lsetp1_{}".format(i))

# %%
square = DomainTypeArray((NEG, NEG, NEG))
with TaskManager():
    square.Compress(lsetsp1)
    boundary = square.Boundary()
    boundary.Compress(lsetsp1)

mlci = MultiLevelsetCutInfo(mesh, lsetsp1)

# %%
# Element and degrees-of-freedom markers
els_if_singe = {dtt: BitArray(mesh.ne) for dtt in boundary}
facets_gp = BitArray(mesh.nedge)

hasneg = mlci.GetElementsWithContribution(square)
Draw(BitArrayCF(hasneg), mesh, "hasneg")

# %% Finite element space
Vhbase = VectorH1(mesh, order=k, dirichlet="left", dgjumps=True)
Vh = Restrict(Vhbase, hasneg)
gfu = GridFunction(Vh)

hasif = mlci.GetElementsWithContribution(boundary)
Draw(BitArrayCF(hasif), mesh, "hasif")

for i, (dtt, els_bnd) in enumerate(els_if_singe.items()):
    els_bnd[:] = mlci.GetElementsWithContribution(dtt)
    Draw(BitArrayCF(els_bnd), mesh, "els_if_singe" + str(i))

facets_gp = GetFacetsWithNeighborTypes(mesh, a=hasneg, b=hasif,
                                       use_and=True)

els_gp = GetElementsWithNeighborFacets(mesh, facets_gp)
Draw(BitArrayCF(els_gp), mesh, "gp_elements")

# %% Bilinear and linear forms of the weak formulation
u, v = Vh.TnT()
h = specialcf.mesh_size
normals = square.GetOuterNormals(lsetsp1)

# Set up the integrator symbols
dx = dCut(lsetsp1, square, definedonelements=hasneg)
ds = {dtt: dCut(lsetsp1, dtt, definedonelements=els_if_singe[dtt])
      for dtt in boundary}
dw = dFacetPatch(definedonelements=facets_gp)


def eps:
    return Sym(Grad(u))


E, nu = 210, 0.2
mu = E / 2 / (1+nu)
lam = E * nu / ((1+nu)*(1-2*nu))

# Construct integrator
a = RestrictedBilinearForm(Vh, facet_restriction=facets_gp, check_unused=False)
a += 2 * mu * InnerProduct(eps(u), eps(v)) * dx + lam * div(u) * div(v) * dx
a += mu * gamma_s / (h**2) * (u - u.Other()) * (v - v.Other()) * dw

f = LinearForm(Vh)
f += rhs * v * dx


# Assemble and solve the linear system
f.Assemble()
a.Assemble()

gfu.vec.data = a.mat.Inverse(Vh.FreeDofs()) * f.vec

# %% Draw webgui
ind = IfPos(lsetsp1[0], 1, IfPos(lsetsp1[1], 1, IfPos(lsetsp1[2], 1, -1)))
DrawDC(ind, gfu[1], 0, mesh, "u")

# %% pyvista vis
try:
    import pyvista
    vtk = VTKOutput(mesh, coefs=[*lsetsp1, 33*CF((gfu[0], gfu[1], 0))],
                    names=["lset1", "lset2", "lset3", "u"], filename="vtkout",
                    subdivision=3, floatsize="single", legacy=False).Do()
    pyvista.set_jupyter_backend('client')
    visobj = pyvista.read('vtkout.vtu')
    plot = pyvista.Plotter()
    arrows = visobj.glyph(scale="u", orient="u", factor=0.15, tolerance=0.03)
    clipped0 = visobj.clip_scalar("lset1", invert=True) \
                     .clip_scalar("lset2", invert=True) \
                     .clip_scalar("lset3", invert=True)
    visobj2 = visobj.warp_by_vector(factor=0.15)
    clipped = visobj2.clip_scalar("lset1", invert=True) \
                     .clip_scalar("lset2", invert=True) \
                     .clip_scalar("lset3", invert=True)
    plot.add_mesh(clipped0, label='reference conf.')
    plot.add_mesh(clipped, label='deformed conf.')
    plot.add_mesh(arrows, label='deform. direction')
    plot.add_legend()
    plot.show()
except ImportError:
    print("pyvista not installed")
    pass
