# Example:
# L2-approximation of a piecewise smooth function which is 
# discontinuous across an implicitly prescribed interface.
#
# The interface is prescribed with the coef "lset",
# which prescribes a simple circle with radius R (=0.4)
#
# Solves the L2 - problem:
#
#   \sum_{i=1,2} \int_{\Omega_i} u_h · v_h dx
# = \sum_{i=1,2} \int_{\Omega_i} u   · v_h dx
#
# Details:
#
# → domains \Omega_1 and \Omega_2: described by coef "lset"
# → FE Space: XH1FESpace, i.e. 
#             standard h1fespace + XFE enrichment
#             use standard h1fespace-flags to adjust space
#             an additional flag "empty" can be used to 
#             avoid the enrichment while still keeping the 
#             composite integration rules. 
#
# Things to try here:
#   1. comment in the shapetester line and observe the shape 
#      functions. To see all shape functions increase 
#      subdivisions to at least 3. Also disable automatic 
#      scaling. Instead use minval=0 and maxval=1.
#      Afterwards comment out the shapetester line.
#   2. refine the mesh and observe the convergence order (L2)
#   3. set flag "empty" (this removes the enrichted functions)
#      for xh1fespace, refine the mesh and observe the 
#      convergence order (L2). Afterwards remove "empty" again.
#   4. set "order" to 2 and try 2. and 3. again, make sure to 
#      increase "ref_space" in xh1fespace and "intorder" in 
#      xdifference
#   5. set "order" back to 1, set the preconditioner to "local" 
#      and add the "-test"-flags, refine the mesh several times
#      and observe the performance
#   6. set "order" to 2 and try 5. again
#


# load geometry
geometry = d1_approx.in2d

# and mesh
mesh = d1_approx.vol.gz

#load xfem-library
shared = libngsxfem_xfem
shared = libngsxfem_py

pymodule = d1_approx

define constant heapsize = 1e9

define constant R = 0.4
define constant one = 1.0

# interface description as zero-level
define coefficient lset
( sqrt(x*x+y*y) - R),       

# solution u- in 'inner' domain (Omega_1)
define coefficient solneg
0.5,

# solution u+ in 'outer' domain (Omega_2)
define coefficient solpos
(sin((x*x+y*y-R*R))),

# use an "extended" continuous finite element space
# you may change the order here
define fespace fescomp
       -type=xstdfespace
       -type_std=h1ho
       -order=1
#       -ref_space=5
#       -empty

#update "extended" part of XFE space:
numproc informxfem npix 
        -xstdfespace=fescomp
        -coef_levelset=lset


define gridfunction u -fespace=fescomp

# integration on sub domains
define linearform f -fespace=fescomp
xsource solneg solpos

# integration on sub domains
define bilinearform a -fespace=fescomp -symmetric -linearform=f
xmass 1.0 1.0

# for the evaluation of gradients
define bilinearform a_d -fespace=fescomp -nonassemble
laplace one -comp=1

numproc drawflux npdf -solution=u -bilinearform=a_d -label=grad

#define preconditioner c -type=local -bilinearform=a -test #-block
define preconditioner c -type=direct -bilinearform=a -inverse=pardiso #-test
#define preconditioner c -type=bddc -bilinearform=a -test -block
#define preconditioner c -type=multigrid -bilinearform=a -test #-smoother=block

#numproc shapetester npst -gridfunction=u #-comp=0

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6 

# evaluate l2-error (difference between prescribed solution and numerical solution)
numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        -levelset=lset
        -intorder=3

numproc visualization npviz -scalarfunction=u #-comp=0
