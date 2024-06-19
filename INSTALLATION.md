# Installation instructions
## Installation of `ngsxfem`
We provide two (and a half) ways to setup `ngsxfem` on your machine:
* [Building/Installing through pip](#installation-through-pip-pip3)
* [Building/Installing from sources](#installation-from-source)
* [Running a docker image](#docker-container)

Below we discuss these installation steps in more detail. If you observe any problems with the installation, you can contact us through the [github issue tracker](https://github.com/ngsxfem/ngsxfem/issues) or the [NGSolve user forum](https://ngsolve.org/forum/index).

Note that your xfem installation depends on your NGSolve installation. 
If you are already NGSolve user and have NGSolve installed, you should choose your xfem-installation based on your NGSolve installation. 
The recommended options are:
 * If you have a release version of NGSolve installed through pip, you should install the release version of xfem through pip as well. 
 * If you have a pre-release are self-built version of NGSolve, you should install xfem through the source wheel. 


## Installation through `pip` (`pip3`)
Table of contents for `pip` install:
  * [Releases](#releases)
  * [Development version](#development-version)
---

We publish a several distribution of the releases of `ngsxfem` on PyPI for installation via `pip3`. 

### Releases
Installation of `ngsxfem` releases through `pip` is carried with the command
``` 
pip3 install xfem
```
You may add standard `pip` options such as `--user`, `--upgrade` and/or `--verbose` or specify a concrete version, e.g. by replacing `xfem` with `xfem==1.4.2104`.
### Development version

For pre-release versions the recommended version is to build from 
source through a source whell.  

``` 
pip3 install xfem --no-binary --no-build-isolation --pre
```
Alternatively, you can try using the binaries also for prelease versions, but compatibility cannot be guaranteed.

## Installation from source

You can also directly install (through `pip`) the latest version of ngsxfem which is the development version in the master branch on github using
``` 
pip3 install git+https://github.com/ngsxfem/ngsxfem.git@master
```
or clone and call cmake manually. You will require (for the latter approach at least) an installed version of NGSolve. This can be pip-installed or from sources. Both should work.

## Docker container
A convenient and reproducible way to set up `ngsxfem` is the usage of [a docker image](https://hub.docker.com/r/ngsxfem/ngsxfem) that we provide here:
<https://hub.docker.com/r/ngsxfem/ngsxfem>.
Installation of `docker` on the common platforms is described [here](https://docs.docker.com/get-docker/). After installation the `Docker daemon` has to be started. This can either be done on boot or manually. In most Linux distributions the command for the latter is either `systemctl start docker` or `service docker start`.

Assuming an installed `docker` with running docker daemon, you can spawn into the `ngsxfem` image with
``` 
docker run -i -t ngsxfem/ngsxfem:latest /bin/bash
```
### Jupyter through docker
To directly spawn a jupyter server in the docker that you can access from a browser start
``` 
docker run -p 8888:8888 ngsxfem/ngsxfem-jupyter
```
and open a browser on your host machine (not in the docker) and paste in the URL that you obtain in the terminal. All computations are carried out in the docker and passed through your browser for interaction. You will have the jupyter tutorial files from the docker container available to work with. Note that changes will not be persistent in the image. To work on local files (with persistent changes) mount a local directory to the docker container, e.g.
``` 
docker run -p 8888:8888 -v ${PWD}:/home/jovyan ngsxfem/ngsxfem-jupyter
```
