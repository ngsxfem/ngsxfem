FROM ngsxfem/ngsolve:v6.2.2105

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

USER root
RUN apt-get install -y cmake g++
USER ${NB_USER}
        
WORKDIR ${HOME}

## manual build based on local directory:
#RUN mkdir ngsxfem
#COPY . ${HOME}/ngsxfem/
#USER root
#RUN chown -R ${NB_UID} ${HOME}/ngsxfem
#USER ${NB_USER}

#WORKDIR ${HOME}/ngsxfem
#RUN ls -al
#RUN pip3 install . --user --upgrade --verbose

## build based on github release:
RUN git clone -b v2.0.2105 --single-branch https://github.com/ngsxfem/ngsxfem.git ngsxfem
RUN pip3 install git+https://github.com/ngsxfem/ngsxfem.git@v2.0.2105 --user --upgrade --verbose
                
RUN python3 -c "import ngsolve; import xfem"        
                
