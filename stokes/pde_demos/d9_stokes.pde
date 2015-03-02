# and mesh
mesh = d9_stokes.vol.gz

#load xfem-library and python-bindings
shared = libngsxfem_xfem                                       
shared = libngsxfem_xstokes                                       
shared = libngsxfem_levelset
# pymodule = d7_stokes
pymodule = instatstokes

pynumproc npInstatStokes np3

