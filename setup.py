
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

#CMAKE_ARGS = ["-DCMAKE_OSX_DEPLOYMENT_TARGET=10.9",
#        "-DCMAKE_PREFIX_PATH=/Users/mhochsteger/Documents/Netgen.app/Contents/Resources",
#        "-DCMAKE_OSX_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk",
#        "-DBUILD_NGSOLVE=OFF",
#        "-DSETUP_PY=ON",
#        "-DUSE_CCACHE=ON"
#        ]


setup(name="xfem",
      version="1.1.0",
      description="",
      packages=["xfem", "xfem.lsetcurv", "xfem.lset_spacetime", "xfem.utils"],
      package_dir={"xfem" : "python",
                   "xfem.lsetcurv" : "lsetcurving",
                   "xfem.lset_spacetime" : "spacetime",
                   "xfem.utils" : "utils"},
      cmake_args=CMAKE_ARGS)

