
# load geometry
geometry = d9_stokes.in2d                                        
# and mesh
mesh = d9_stokes_finer.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
# pymodule = d7_stokes
pymodule = instatstokes

constant heapsize = 1e9

constant R = 0.006

constant x0 = 0
constant x1 = 0.003

constant y0 = 0

coefficient initial_lset
( sqrt((x-x0)*(x-x0)+4*(y-y0)*(y-y0)) - R),

# define coefficient initial_lset
# (
# (sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))-R*R < sqrt((x-x1)*(x-x1)+(y-y0)*(y-y0))-R*R) * (sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))-R*R) 
# +
# (sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))-R*R > sqrt((x-x1)*(x-x1)+(y-y0)*(y-y0))-R*R) * (sqrt((x-x1)*(x-x1)+(y-y0)*(y-y0))-R*R)
# ),

pynumproc npInstatStokes np3


# define coefficient velocity ( (uvp.1, uvp.2) )
# numproc draw npd1 -coefficient=velocity -label=velocity

# define coefficient pressure ( (uvp.3) )
# numproc draw npd2 -coefficient=pressure -label=pressure

# numproc draw npd3 -coefficient=lset -label=levelset

# numproc visualization npviz 
#         -scalarfunction=levelset
#         # -vectorfunction=velocity
#         -minval=0 -maxval=0 
#         -nolineartexture -deformationscale=1 -subdivision=3

