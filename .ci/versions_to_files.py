import os
import sys

from subprocess import check_output
def path_to_version(path):
    git_version = check_output(['git', 'describe', '--tags'], cwd=path).decode('utf-8').strip()
    version = git_version[1:].split('-')
    if len(version)>2:
        version = version[:2]
    if len(version)>1:
        version = '.post'.join(version) + '.dev'
    else:
        version = version[0]
    return version

version = path_to_version(".")
ngsolve_version = path_to_version("external_dependencies/ngsolve")

f = open("ngsxfem.version", "w")
f.write(version)
f.close()
f = open("ngsolve.version", "w")
f.write(ngsolve_version)
f.close()
