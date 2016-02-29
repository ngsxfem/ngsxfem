
geometry = cube_simple.geo
#mesh = cube.vol.gz
mesh = cube_simple.vol.gz


#load xfem-library and python-bindings
shared = libngsxfem_xfem
shared = libngsxfem_tracefem
shared = libngsxfem_lsetcurving


constant heapsize = 1e9

constant R = 1.0
constant one = 1.0
        
####levelsets
coefficient lset
( sqrt(x*x+y*y+z*z) - R),

fespace fes_lset_ho -type=h1ho -order=4
gridfunction lset_ho -fespace=fes_lset_ho
numproc setvalues npsv1 -gridfunction=lset_ho -coefficient=lset
        
###########################################################################
########################### quasi-normal field ############################
###########################################################################
fespace fes_normal -type=h1ho -order=4 -vec #-dirichlet=[1,2,3,4,5,6]
gridfunction qn -fespace=fes_normal
numproc setvalues npsv3 -gridfunction=qn -coefficient=grad_lset_ho

        
###########################################################################
########################### P1 ############################################
###########################################################################
        
### project this level set function into a finite element space of order 1
fespace fes_lset_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_lset_p1
numproc interpolatep1 npipp1b -gridfunction_ho=lset_ho -gridfunction_p1=lset_p1

### determine the deformation 
fespace fes_deform -type=h1ho -order=4 -vec -dirichlet=[1,2,3,4,5,6]
gridfunction deform -fespace=fes_deform

numproc projectshift nppsh -levelset=lset_ho -levelset_p1=lset_p1 -deform=deform -quasinormal=qn -lset_lower_bound=0.0 -lset_upper_bound=0.0 -threshold=50.5
                                                

###########################################################################
########################### deformation ###################################
###########################################################################

# numproc calcerrors npcalcerr -levelset_ho=lset_ho -levelset_p1=lset_p1 -deform=deform
#                 -lset_lower_bound=0
#                 -lset_upper_bound=0

numproc setdeformation npudef -gridfunction=deform
                
define fespace fesh1
       -type=h1ho
       -dirichlet=[1,2,3,4,5,6]
       # -dgjumps
       -order=4

define fespace tracefes
       -type=xfespace
       -type_std=h1ho
       -ref_space=0
       -dirichlet=[1,2,3,4,5,6]
       -trace
       # -dgjumps
       -order=4

numproc informxfem npix
        -xfespace=tracefes
        -fespace=fesh1
        -coef_levelset=lset_p1

gridfunction gf_u -fespace=tracefes

coefficient normal
(x,y,z),

bilinearform a -fespace=tracefes
tracemass 1.0
tracelaplacebeltrami 1.0
#tracelaplace 1.0
# lo_traceghostpenalty 0.01
#sec_traceghostpenalty 0.001

#normallaplacetrace 0.1 normal
        
linearform f -fespace=tracefes
# tracesource (sin(pi*z)*(1+pi*pi*(1-z*z*z))+cos(pi*z)*4*pi*z)
# tracesource (z*z+6*z*z-2)#solution u=z**2
#tracesource (z) #solution u=z ???
tracesource (sin(pi*z)*(pi*pi*(1-z*z)+1)+cos(pi*z)*2*pi*z) #solution u=sin(pi*z)
# tracesource (sin(pi*z))

# coefficient u_sol
# (z*z),

coefficient u_sol
(sin(pi*z)),

# coefficient gradu_gamma_sol
# ((-x*pi*z*cos(pi*z)),(-y*pi*z*cos(pi*z)),(pi*cos(pi*z)-z*pi*z*cos(pi*z))),
        
define preconditioner c -type=local -bilinearform=a -test #-block
#efine preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test
# define preconditioner c -type=bddc -bilinearform=a -test #-block

numproc bvp npbvp -gridfunction=gf_u -bilinearform=a -linearform=f -solver=cg -preconditioner=c -maxsteps=100000 -prec=1e-7
# define preconditioner c -type=direct -bilinearform=a -inverse=pardiso -test


bilinearform evalu -fespace=tracefes -nonassemble
exttrace 1.0

numproc drawflux npdf -solution=gf_u -bilinearform=evalu -applyd -label=u

numproc visualization npviz -scalarfunction=u -vectorfunction=deform
    -nolineartexture -deformationscale=1.0 -subdivision=2

       
numproc tracediff3d npxd 
        -gridfunction=gf_u 
        -coef=u_sol
        -intorder=6

                
numproc unsetdeformation npudef
        
# numproc vtkoutput npout -filename=simple
#         -coefficients=[lset,lset_ho,lset_p1,deform,gf_u,u_sol]
#         -fieldnames=[lset,lsetho,lsetp1,deformation,u,u_sol]
#         -subdivision=2

numproc levelsetrefine nplsref -levelset=lset_p1

# numproc markinterface npmi -fespace=tracefes
        