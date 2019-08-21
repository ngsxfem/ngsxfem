
try:
    import skbuild
except ImportError:
    print("")
    print("***********************************************************************************")
    print("")
    print("skbuild not found, please execute:")
    print("pip3 install scikit-build")
    print("")
    print("***********************************************************************************")
    print("")

from skbuild import setup

CMAKE_ARGS = [ "-DBUILD_NGSOLVE=OFF" ]

# CMAKE_ARGS.append( "-DCMAKE_OSX_DEPLOYMENT_TARGET=10.9" )
# CMAKE_ARGS.append( "-DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk" )

setup(name="xfem",
      version="1.1.0",
      description="",
      packages=["xfem", "xfem.lsetcurv", "xfem.lset_spacetime", "xfem.utils"],
      package_dir={"xfem" : "python",
                   "xfem.lsetcurv" : "lsetcurving",
                   "xfem.lset_spacetime" : "spacetime",
                   "xfem.utils" : "utils"},
      cmake_args=CMAKE_ARGS)

