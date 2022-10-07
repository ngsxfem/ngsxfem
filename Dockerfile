FROM ngsxfem/ngsolve:latest

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

USER root
RUN apt-get install -y cmake g++
USER ${NB_USER}
        
WORKDIR ${HOME}

RUN git clone -b v2.0.2204 --single-branch https://github.com/ngsxfem/ngsxfem.git ngsxfem
RUN pip3 install git+https://github.com/ngsxfem/ngsxfem.git@v2.0.2204 --user --upgrade --verbose
                
RUN python3 -c "import ngsolve; import xfem"        
                
