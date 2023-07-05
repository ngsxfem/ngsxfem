#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import platform
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion
from distutils.sysconfig import get_python_lib


from subprocess import check_output
def path_to_version(path):
    suffix = None
    git_version = check_output(['git', 'describe', '--tags'], cwd=path).decode('utf-8').strip()
    version = git_version[1:].split('-')
    main_version = version[0].split('.')
    if (len(main_version)>3):
        suffix = main_version[3]
        main_version = main_version[0:3]
        version[0] = main_version[0]+"."+main_version[1]+"."+main_version[2]
    if len(version)>2:
        version = version[:2]
    if len(version)>1:
        version = '.post'.join(version)
        if suffix:
            version += "." + suffix
    else:
        version = version[0]
    return version

try:
    version = path_to_version(".")
except:
    version = "2.1.2302.dev1"

try:
    ngsolve_version = path_to_version("external_dependencies/ngsolve")
except:
    ngsolve_version = None

print("ngsxfem_version =", version)
print("ngsolve_version =", ngsolve_version)
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        
class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))
        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")
        for ext in self.extensions:
            self.build_extension(ext)
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                      '-DCMAKE_CXX_COMPILER=ngscxx',
                      '-DCMAKE_LINKER=ngsld',
                      '-DBUILD_STUB_FILES=ON',
                      '-DCHECK_NGSOLVE_VERSION=OFF', # temporary
                      '-DBUILD_NGSOLVE=OFF']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]
        if platform.system() == "Windows":
            #not expected to work... (but who knows..)
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']


           
        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())

        if 'PYDIR' in os.environ:
            cmake_args += [f'-DCMAKE_PREFIX_PATH={os.environ["PYDIR"]}/..']
            cmake_args += [f'-DPYTHON_EXECUTABLE={os.environ["PYDIR"]}/python3']
            cmake_args += [f'-DPYTHON_LIBRARY={os.environ["PYDIR"]}/../lib']
            cmake_args += [f'-DPYTHON_INCLUDE_DIR={os.environ["PYDIR"]}/../include']

        
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
        
        subprocess.check_call(['mv', 'ngsxfem_py.so', 'xfem'], cwd=self.build_lib)

name = "xfem"
ngsolve_name = "ngsolve"
#if 'NETGEN_ARCH' in os.environ and os.environ["NETGEN_ARCH"] == "avx2":
#    name = "xfem-avx2"
#    ngsolve_name = "ngsolve-avx2"

if ngsolve_version:    
    print("require", ngsolve_name+">="+ngsolve_version)
    install_requires = [ngsolve_name+">="+ngsolve_version]
else:
    install_requires = []

setup(
    name=name,
    version=version, #'2.0.dev2',
    author='Christoph Lehrenfeld',
    author_email='lehrenfeld@math.uni-goettingen.de',
    description='(ngs)xfem is an Add-on library to Netgen/NGSolve for unfitted/cut FEM.',
    long_description='(ngs)xfem is an Add-on library to Netgen/NGSolve which enables the use of unfitted finite element technologies known as XFEM, CutFEM, TraceFEM, Finite Cell, ... . ngsxfem is an academic software. Its primary intention is to facilitate the development and validation of new numerical methods.',
    url="https://github.com/ngsxfem/ngsxfem",
    install_requires=install_requires,
    ext_modules=[CMakeExtension('ngsxfem_py')],
    cmdclass=dict(build_ext=CMakeBuild),
    packages=["xfem"],
    package_dir={"xfem": "python"},
    python_requires='>=3.8',
)
