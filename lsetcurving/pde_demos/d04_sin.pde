geometry = d01_testgeom.in2d
mesh = d01_01.vol.gz
#mesh = d04_refined.vol.gz
#mesh = d04_unstr.vol.gz
# mesh = d99_testgeom_unstr.vol.gz
shared = libngsxfem_xfem
shared = libngsxfem_lsetcurving

flags tracer = -max_size=0
        
define constant order_deform = "5"
define constant order_qn = "1"
define constant order_lset = "5"

constant levelset_lower_bound = "-0.0"
constant levelset_upper_bound = "0.0"
constant lset_lower_bound = -0
constant lset_upper_bound = 0

        
constant one = 1.0
constant threshold = "0.1"
        #0.5/$(order_deform)
# constant threshold = 0.1
        
define constant heapsize = 2e9
constant R = 2.0/3.0

define constant omega = 8.0
define constant omega2 = (1.2*pi) 
        
# # (y-0.001*sin(8*pi*(x))),
# define coefficient lset
# (
#  sqrt(y*y+x*x)-(R+0.05*sin(omega*atan2(x,y))
# )
        
# (sqrt(x*x+y*y)-0.1),

# define coefficient lset
# ((sqrt((x)*(x)+(y)*(y))-0.6666666666666666666666)*(1+0.5*sin(8*pi*x)*cos(8*pi*y))),



constant x0 = 0.0
constant y0 = 0.0
constant R0 = 0.5

constant x1 = -0.25
constant y1 = 0.25
constant R1 = 0.5

constant x2 = 0.5
constant y2 = -0.5
constant R2 = 0.25

#geometry constants


define coefficient lset
(
 sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))-(R0+0.1*sin(omega*atan2(x-x0,y-y0)))
)
        
define coefficient obj1
(
 sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))-(R1+0.1*sin(omega*atan2(x-x1,y-y1)))
)
                        
define coefficient r1
(sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))),

# define coefficient obj1
# (r1-R1),

define coefficient r2
(sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2))),
        
define coefficient obj2
(r2-R2),

# define coefficient lset
# (
#           ((obj2*obj2) > (obj1*obj1)) * obj1 
#         + ((obj1*obj1) > (obj2*obj2)) * obj2
# )

        
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

# ###########################################################################
# ########################### quasi-normal field ############################
# ###########################################################################
# fespace fes_normal -type=h1ho -order=$(order_qn) -vec #-dirichlet=[1,2,3,4,5,6]
# gridfunction qn -fespace=fes_normal
        
# bilinearform mqn -fespace=fes_normal
# mass 1.0 --comp=1        
# mass 1.0 --comp=2
# #mass 1.0 --comp=1        

# coefficient test1
# (grad_lset_ho*(1,0))

# coefficient test2
# (grad_lset_ho*(0,1))

# # coefficient test1
# # ((x+1))

# # coefficient test2
# # ((y+1))
                        
# linearform fqn -fespace=fes_normal
# source test1 --comp=1
# source test2 --comp=2
# # source (grad_lset_ho*(0,1)) --comp=2


# define preconditioner cqn -type=local -block -blocktype=7 -bilinearform=mqn
              
# numproc bvp npbvpqn -gridfunction=qn -bilinearform=mqn -linearform=fqn -solver=cg -preconditioner=cqn -maxsteps=1000 -prec=1e-6 # -print
        
###########################################################################
########################### quasi-normal field ############################
###########################################################################
        
### project this level set function into a finite element space of order 1
fespace fes_p1 -type=h1ho -order=1
gridfunction lset_p1 -fespace=fes_p1

numproc interpolatep1 npipp1b -gridfunction_ho=lset_ho -gridfunction_p1=lset_p1

### determine the deformation 
fespace fes_deform -type=h1ho -order=$(order_deform) -vec -dirichlet=[1,2,3,4]
        
gridfunction deform -fespace=fes_deform
                                        
# bilinearform mdeform -fespace=fes_deform -symmetric
# restrictedmass_4 one lset_p1 lset_lower_bound lset_upper_bound --comp=1        
# restrictedmass_4 one lset_p1 lset_lower_bound lset_upper_bound --comp=2
# #mass 1.0 --comp=1        

# linearform fdeform -fespace=fes_deform
# #shiftsource_5 lset_p1 lset_ho threshold lset_lower_bound lset_upper_bound
# shiftsource_6 lset_p1 lset threshold lset_lower_bound lset_upper_bound qn

# define preconditioner cdeform -type=direct -bilinearform=mdeform
              
# numproc bvp npbvpdeform -gridfunction=deform -bilinearform=mdeform -linearform=fdeform -solver=cg -preconditioner=cdeform -maxsteps=1000 -prec=1e-6 # -print

numproc projectshift nppsh -levelset=lset_ho -levelset_p1=lset_p1 -deform=deform -quasinormal=qn -lset_lower_bound=$(lset_lower_bound) -lset_upper_bound=$(lset_upper_bound) -threshold=0.1
        #$(threshold)
        
numproc unsetdeformation npudef

numproc levelsetrefine nplsref -levelset=lset_p1

numproc calcerrors npcalcerr -levelset_ho=lset -levelset_p1=lset_p1 -quasinormal=qn -deform=deform
                -lset_lower_bound=$(levelset_lower_bound)
                -lset_upper_bound=$(levelset_upper_bound)


numproc visualization npvis -scalarfunction=lset_p1 -vectorfunction=deform -deformationscale=1 -subdivision=0  -minval=0.0 -maxval=0.0 -nolineartexture

# numproc vtkoutput npout -filename=d04_output_sd5_
#         -coefficients=[lset]
#         -fieldnames=[levelset]
#         -subdivision=5
        
#        -gridfunctions=[gf_negpos.1,gf_negpos.2]

        