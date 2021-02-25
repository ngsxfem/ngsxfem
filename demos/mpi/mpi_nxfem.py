"""
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
import os
from mpi4py import MPI

import netgen.meshing
from netgen.geom2d import SplineGeometry
from ngsolve import *
from ngsolve import ngs2petsc
from xfem import *
from xfem.lsetcurv import *

from math import pi


comm = MPI.COMM_WORLD
rank = comm.rank
np = comm.size

ngsglobals.msg_level = 10


# -------------------------------- PARAMETERS ---------------------------------
# Domain corners
ll, ur = (-1.5, -1.5), (1.5, 1.5)
# Mesh diameter
maxh = 0.2
# Finite element space polynomial order
order = 2
# diffusion coefficients for the sub-domains (NEG/POS):
alpha = [1.0, 2.0]
# Nitsche penalty parameter
lambda_nitsche = 20
# Write VTK output for visualisation
do_vtk = True


# ----------------------------------- MAIN ------------------------------------
# Geometry and Mesh
geo = SplineGeometry()
geo.AddRectangle(ll, ur, bc=1)

if rank == 0:
    ngmesh = geo.GenerateMesh(maxh=maxh)
    if np > 1:
        ngmesh.Distribute(comm)
else:
    ngmesh = netgen.meshing.Mesh.Receive(comm)
    ngmesh.SetGeometry(geo)

mesh = Mesh(ngmesh)

# Level set function of the domain and higher-order machinery
levelset = (sqrt(sqrt(x**4 + y**4)) - 1)

lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000,
                                      discontinuous_qn=True)
deformation = lsetmeshadap.CalcDeformation(levelset)
lsetp1 = lsetmeshadap.lset_p1
Draw(lsetp1, mesh, "lsetp1")

mesh.SetDeformation(deformation)

# Gather information of cut elements
ci = CutInfo(mesh, lsetp1)

# Extended finite element space
Vh = H1(mesh, order=order, dirichlet=[1, 2, 3, 4])
Vhx = XFESpace(Vh, ci)

VhG = FESpace([Vh, Vhx])
print("unknowns in extended FESpace for rank", rank, ":", VhG.ndof)

# manufactured solution and corresponding r.h.s. data CoefficientFunctions:
# Manufactured solution and corresponding r.h.s. data CoefficientFunctions:
r44 = (x**4 + y**4)
r41 = sqrt(sqrt(x**4 + y**4))
r4m3 = (1.0 / (r41 * r41 * r41))
r66 = (x**6 + y**6)
r63 = sqrt(r66)
r22 = (x**2 + y**2)
r21 = sqrt(r22)
solution = [1 + pi / 2 - sqrt(2.0) * cos(pi / 4 * r44), pi / 2 * r41]
coef_f = [(-1 * sqrt(2) * pi * (pi * cos(pi / 4 * (r44)) * (r66)
                                + 3 * sin(pi / 4 * (r44)) * (r22))),
          (-2 * pi * 3 / 2 * (r4m3) * (-(r66) / (r44) + (r22)))]

# Coefficients:
n = 1.0 / grad(lsetp1).Norm() * grad(lsetp1)
h = specialcf.mesh_size

# the cut ratio extracted from the cutinfo-class
kappa = (CutRatioGF(ci), 1.0 - CutRatioGF(ci))
# Nitsche stabilization parameter:
stab = lambda_nitsche * (alpha[1] + alpha[0]) / h

# expressions of test and trial functions:
u_std, u_x = VhG.TrialFunction()
v_std, v_x = VhG.TestFunction()

u = [u_std + op(u_x) for op in [neg, pos]]
v = [v_std + op(v_x) for op in [neg, pos]]

gradu = [grad(u_std) + op(u_x) for op in [neg_grad, pos_grad]]
gradv = [grad(v_std) + op(v_x) for op in [neg_grad, pos_grad]]

average_flux_u = sum([- kappa[i] * alpha[i] * gradu[i] * n for i in [0, 1]])
average_flux_v = sum([- kappa[i] * alpha[i] * gradv[i] * n for i in [0, 1]])

# Integration domains for integration on negative/positive sub-domains
# and on the interface: Here, the integration is (geometrically) exact
# if the "levelset"-argument is a piecewise (multi-)linear function. The
# integration order is chosen according to the arguments in the
# multi-linear forms (but can be overwritten with "force_intorder" in
# the integration domain). If the "levelset"-argument is not a
# (multi-)linear function, you can use the "subdivlvl" argument to add
# additional refinement levels for the geometry approximation.
lset_neg = {"levelset": lsetp1, "domain_type": NEG, "subdivlvl": 0}
lset_pos = {"levelset": lsetp1, "domain_type": POS, "subdivlvl": 0}
lset_if = {"levelset": lsetp1, "domain_type": IF, "subdivlvl": 0}

# bilinear forms:
a = BilinearForm(VhG, symmetric=True, nonsym_storage=True)
# l.h.s. domain integrals:
a += SymbolicBFI(levelset_domain=lset_neg, form=alpha[0] * gradu[0] * gradv[0])
a += SymbolicBFI(levelset_domain=lset_pos, form=alpha[1] * gradu[1] * gradv[1])
# Nitsche integrals:
a += SymbolicBFI(levelset_domain=lset_if, form=average_flux_u * (v[0] - v[1])
                 + average_flux_v * (u[0] - u[1])
                 + stab * (u[0] - u[1]) * (v[0] - v[1]))

f = LinearForm(VhG)
# r.h.s. domain integrals:
f += SymbolicLFI(levelset_domain=lset_neg, form=coef_f[0] * v[0])
f += SymbolicLFI(levelset_domain=lset_pos, form=coef_f[1] * v[1])

# Solution vector
gfu = GridFunction(VhG)
gfu.components[0].Set(solution[1], BND)

# # setting up matrix and vector
# c = Preconditioner(a, 'hypre') # very good for low order but fine meshes
# c = Preconditioner(a, 'bddc') # very good for high order but moderate meshes
# c = Preconditioner(a, 'direct', inverse="masterinverse") # only for small run
# c = Preconditioner(a, 'gamg') # petsc gamg
if np == 1:
    c = Preconditioner(a, 'bddc', coarsetype="h1amg")  # only for serial runs
else:
    c = Preconditioner(a, 'bddc', coarsetype="gamg")

a.Assemble()
f.Assemble()

rhs = gfu.vec.CreateVector()
rhs.data = f.vec - a.mat * gfu.vec
update = gfu.vec.CreateVector()
update.data = solvers.CG(mat=a.mat, pre=c.mat, rhs=rhs,
                         tol=1e-6, maxsteps=100, printrates=comm.rank == 0)
gfu.vec.data += update

uh = [gfu.components[0] + op(gfu.components[1]) for op in [neg, pos]]
err_sqr_coefs = [(uh[i] - solution[i]) * (uh[i] - solution[i]) for i in [0, 1]]
# err_sqr_coefs = [uh[i]**2 for i in [0,1]]

l2error = sqrt(Integrate(levelset_domain=lset_neg, cf=err_sqr_coefs[0],
                         mesh=mesh, order=2)
               + Integrate(levelset_domain=lset_pos, cf=err_sqr_coefs[1],
                           mesh=mesh, order=2))

# l2error = Integrate(cf=IfPos(lsetp1,err_sqr_coefs[1],err_sqr_coefs[0]),
#                     mesh=mesh, order=2)

# print(Integrate(levelset_domain = lset_neg, cf=1, mesh=mesh, order=2))
# print(Integrate(levelset_domain = lset_pos, cf=1, mesh=mesh, order=2))
# print(Integrate(levelset_domain = lset_pos, cf=1, mesh=mesh, order=2)
#       + Integrate(levelset_domain = lset_neg, cf=1, mesh=mesh, order=2))
if rank == 0:
    print("L2 error : ", l2error)

# Do VTK output for visualisation
output_path = os.path.dirname(os.path.realpath(__file__)) + "/output"
if rank == 0 and not os.path.exists(output_path):
    os.mkdir(output_path)
comm.Barrier()  # wait until master has created the directory!!

u = [gfu.components[0] + op(gfu.components[1]) for op in [neg, pos]]

if do_vtk:
    vtk_output = [u[0] - solution[0], u[1] - solution[1], lsetp1, deformation]
    vtk_names = ["u1", "u2", "lset", "deformation"]

    vtk = VTKOutput(ma=mesh, coefs=vtk_output, names=vtk_names,
                    filename=output_path + "/vtkout_p" + str(rank) + "_n0",
                    subdivision=1)
    vtk.Do()
comm.Barrier()
