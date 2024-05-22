##################################################
# This file generates a docker image with 
#  xfem installed from source on top of pip-NGSolve.
# Usage:
#    `docker build -t xfem .`
#  in the directory containing the Dockerfile 
#  generates the image.
##################################################

# this file is partially copied from gitlab.gwdg.de/ngsuite/docker-ngsolve
FROM ubuntu:24.04 as prebuild

ARG NB_USER=jovyan
ARG NB_UID=1001

ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN useradd --shell /usr/sbin/nologin \
--uid ${NB_UID} \
${NB_USER}

WORKDIR ${HOME}
RUN mkdir ${HOME}/ngsxfem
RUN mkdir ${HOME}/build

RUN export DEBIAN_FRONTEND=noninteractive
RUN apt-get update

RUN apt-get install -y python3.12-venv python3.12-dev cmake
RUN apt-get update
RUN ln -s /usr/bin/python3.12 /usr/bin/python3

ARG NGS_VERSION=6.2.2402
# ARG XFEM_VERSION=v2.1.2303

# install pip-packages in a virtual environment and set PYTHONPATH to the according venv-directory
#  this avoids "externally managed environment" error
RUN /bin/bash -c "python3 -m venv --clear .venv &&\
                    source ./.venv/bin/activate &&\
                    pip3 install --no-cache-dir ngsolve==${NGS_VERSION} &&\
                    pip3 install --no-cache-dir notebook --no-cache-dir jupyterlab --no-cache-dir numpy --no-cache-dir scipy --no-cache-dir matplotlib --no-cache-dir ipywidgets --no-cache-dir psutil --no-cache-dir pytest --no-cache-dir webgui_jupyter_widgets &&\
                    pip3 list"

RUN apt-get install -y git
RUN apt-get install -y g++

# gitlab/xfem is currently not public
RUN git clone --single-branch https://github.com/ngsxfem/ngsxfem.git

# github/xfem currently needs explicit cmake cxx flags 
RUN cmake -DCMAKE_INSTALL_PREFIX=${HOME}/.venv -DBUILD_NGSOLVE=OFF -DNGSolve_DIR=/home/.venv/lib/cmake/ngsolve -DCMAKE_CXX_FLAGS='-fabi-version=14 -D_GLIBCXX_USE_CXX11_ABI=0' -DCHECK_NGSOLVE_VERSION=OFF ${HOME}/ngsxfem
RUN make install


FROM ubuntu:24.04 as postbuild

ARG NB_USER=jovyan
ARG NB_UID=1001

ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN useradd --shell /usr/sbin/nologin \
    --uid ${NB_UID} \
    ${NB_USER}

WORKDIR ${HOME}
RUN export DEBIAN_FRONTEND=noninteractive
RUN apt-get update

RUN apt-get install -y cmake python3-pip git g++
RUN apt-get update
# RUN ln -s /usr/bin/python3.12 /usr/bin/python3

COPY --from=prebuild ${HOME}/.venv /home/.venv
ENV PYTHONPATH "${PYTHONPATH}:/home/.venv/lib/python3.12/site-packages"
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/home/.venv/lib"

RUN chown -R ${NB_UID} ${HOME}

WORKDIR ${HOME}
RUN python3 -c "import ngsolve; import xfem"

ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]

USER ${NB_USER}
WORKDIR /home/${NB_USER}

CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root" ]