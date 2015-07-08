geometry = d99_testgeom.in2d
mesh = d99_testgeom_unstr.vol.gz
# mesh = d99_testgeom_unstr.vol.gz
shared = libngsxfem_xfem

define constant heapsize = 1e9
constant R = 0.5-1e-12

define coefficient lset
( sqrt(x*x+y*y) - R),


define constant aneg = 2.0
define constant apos = 1.0


define coefficient rhsneg
(8),

define coefficient rhspos
(2/(sqrt(x*x+y*y))),

define coefficient solpos
(1.0-2.0*sqrt(x*x+y*y)),

define coefficient solneg
(1.0/4.0-(x*x+y*y)),


define fespace fescomp
       -type=xstdfespace
       -type_std=h1ho 
       -order=2
       -dirichlet=[1,2]
       -ref_space=0
#       -dgjumps

numproc informxfem npix 
        -xstdfespace=fescomp
        -coef_levelset=lset

define gridfunction u -fespace=fescomp

numproc setvaluesx npsvx -gridfunction=u -coefficient_neg=solneg -coefficient_pos=solpos -boundary #-print
        
######### CURV IT #########
        
fespace fes_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_p1

fespace fes_p2 -type=h1ho -order=2
gridfunction lset_p2 -fespace=fes_p2
                
fespace fes_deform -type=h1ho -order=2 -dim=2
gridfunction deform -fespace=fes_deform
                
numproc xgeomtest npxd 
        -gf_levelset_p1=lset_p1
        -gf_levelset_p2=lset_p2
        -levelset=lset
        -gf_levelset=lset_p2
        -deformation=deform
        # -nocutoff
        -threshold=0.8

######### CURV IT #########

define linearform f -fespace=fescomp # -print
xsource rhsneg rhspos

define bilinearform a -fespace=fescomp -printelmat #-eliminate_internal -keep_internal -symmetric -linea
xlaplace aneg apos
xnitsche_heaviside aneg apos 1.0 1.0 10.0
        
define preconditioner c -type=direct -bilinearform=a

numproc bvp npbvp -gridfunction=u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=1000 -prec=1e-6 # -print

                
numproc draw npdr -coefficient=lset -label=levelset


numproc visualization npvis -scalarfunction=lset_p1 -vectorfunction=deform -deformationscale=1 -subdivision=0 -minval=0 -maxval=0

coefficient zero
0,

numproc xdifference npxd 
        -solution=u 
        -solution_n=solneg
        -solution_p=solpos
        -jumprhs=zero
        -levelset=lset
        -interorder=5
        -henryweight_n=1
        -henryweight_p=1
        -diffusion_n=2
        -diffusion_p=1

coefficient err
(
  ((lset) > 0) * (abs((u)-solpos)) 
 +((lset) < 0) * (abs((u)-solneg)) 
)

numproc draw npdraw -coefficient=err -label=error
        
numproc unsetdeformation npunset

        