FROM ngsxfem/ngsolve:v6.2.2303

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

USER root
RUN apt-get install -y cmake g++
USER ${NB_USER}
        
WORKDIR ${HOME}

ENV LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/usr/local/lib/:/home/jovyan/.local/lib/"
ENV PATH "${PATH}:/home/jovyan/.local/bin"

RUN git clone -b v2.1.2303 --single-branch https://github.com/ngsxfem/ngsxfem.git ngsxfem
RUN pip3 install git+https://github.com/ngsxfem/ngsxfem.git@v2.1.2303 --user --upgrade --verbose --no-dependencies
                
RUN python3 -c "import ngsolve; import xfem"        

CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root" ]

