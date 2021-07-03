FROM ngsxfem/ngsolve:latest

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

USER ${NB_USER}
        
WORKDIR ${HOME}

RUN git clone -b master --single-branch https://github.com/ngsxfem/ngsxfem.git ngsxfem
RUN pip3 install git+https://github.com/ngsxfem/ngsxfem.git@master --user --upgrade --verbose
                
RUN python3 -c "import ngsolve; import xfem"        
                
