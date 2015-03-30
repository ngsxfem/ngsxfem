
# load geometry
geometry = shapetest.in2d                                        
# and mesh
mesh = shapetest.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_py                                       

define constant heapsize = 1e9

define constant R = 0.4s
define constant one = 1.0

# interface description as zero-level
define coefficient lset
( sqrt(x*x+y*y) - R),

define fespace fesh1                                            
       -type=h1ho                                          
       -order=1                                                  

# use an "extended" continuous finite element space
# you may change the order here
define fespace tracefes
       -type=xfespace                                      
       -type_std=h1ho                                          
       -order=1                                                  
       -ref_space=1                                           
#       -empty

#update "extended" part of XFE space:
numproc informxfem npix                                     
        -xfespace=tracefes
        -fespace=fesh1
        -coef_levelset=lset

gridfunction u -fespace=tracefes

numproc shapetester npst -gridfunction=u

numproc visualization npviz -scalarfunction=u 
    -minval=0 -maxval=1 
    -nolineartexture -deformationscale=1 -subdivision=4

#TODO visualize:
# * use xfem visualizer
# or
# * write your own visualization routine