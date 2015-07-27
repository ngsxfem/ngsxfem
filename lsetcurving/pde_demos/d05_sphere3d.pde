geometry = cube.geo
mesh = d05_cube.vol.gz
#mesh = d04_refined.vol.gz
#mesh = d04_unstr.vol.gz
# mesh = d99_testgeom_unstr.vol.gz
shared = libngsxfem_xfem
shared = libngsxfem_lsetcurving

define constant order_deform = "3"
define constant order_qn = "3"
define constant order_lset = "3"

constant levelset_lower_bound = "-0.0"
constant levelset_upper_bound = "0.0"
constant lset_lower_bound = -0
constant lset_upper_bound = 0

        
constant one = 1.0
constant threshold = 10.2 #10.6/$(order_deform)
# constant threshold = 0.1
        
define constant heapsize = 2e9
constant R = 2.0/3.0

# # (y-0.001*sin(8*pi*(x))),
# # (sqrt((x+1)*(x+1)+(y+1)*(y+1))-0.6666666666666666666666),

define coefficient lset
(sqrt((x)*(x)+(y)*(y)+z*z)-0.6666666666666666666666),


#geometry constants
define constant r0 = 0.5+1e-12
define constant omega = 16.0
define constant omega2 = (1.2*pi) 

                
### project level set function into a finite element space of order k
fespace fes_ho -type=h1ho -order=$(order_lset)
gridfunction lset_ho -fespace=fes_ho

numproc setvalues npsv -gridfunction=lset_ho -coefficient=lset

###########################################################################
########################### quasi-normal field ############################
###########################################################################

fespace fes_normal -type=h1ho -order=$(order_qn) -vec #-dirichlet=[1,2,3,4,5,6]
gridfunction qn -fespace=fes_normal
        
numproc setvalues npsv -gridfunction=qn -coefficient=grad_lset_ho 

###########################################################################
########################### quasi-normal field ############################
###########################################################################
        
### project this level set function into a finite element space of order 1
fespace fes_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_p1

numproc interpolatep1 npipp1b -gridfunction_ho=lset_ho -gridfunction_p1=lset_p1

### determine the deformation 
fespace fes_deform -type=h1ho -order=$(order_deform) -vec
        
gridfunction deform -fespace=fes_deform
                                        
bilinearform mdeform -fespace=fes_deform -symmetric
restrictedmass_4 one lset_p1 lset_lower_bound lset_upper_bound --comp=1        
restrictedmass_4 one lset_p1 lset_lower_bound lset_upper_bound --comp=2
restrictedmass_4 one lset_p1 lset_lower_bound lset_upper_bound --comp=3
#mass 1.0 --comp=1        

linearform fdeform -fespace=fes_deform
#shiftsource_5 lset_p1 lset_ho threshold lset_lower_bound lset_upper_bound
shiftsource_6 lset_p1 lset threshold lset_lower_bound lset_upper_bound qn

define preconditioner cdeform -type=direct -bilinearform=mdeform #-inverse=sparsecholesky
              
numproc bvp npbvpdeform -gridfunction=deform -bilinearform=mdeform -linearform=fdeform -solver=cg -preconditioner=cdeform -maxsteps=1000 -prec=1e-6 # -print
                       
numproc unsetdeformation npudef

numproc levelsetrefine nplsref -levelset=lset_p1

numproc calcerrors npcalcerr -levelset_ho=lset -levelset_p1=lset_p1 -quasinormal=qn -deform=deform
                -lset_lower_bound=$(levelset_lower_bound)
                -lset_upper_bound=$(levelset_upper_bound)


numproc visualization npvis -scalarfunction=lset_p1 -vectorfunction=deform -deformationscale=1 -subdivision=3  -minval=0.0 -maxval=0.0 -nolineartexture
        