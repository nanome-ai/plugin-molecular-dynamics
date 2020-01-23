FROM nvidia/cuda AS cuda
# FROM nvidia/opencl as opencl
FROM continuumio/miniconda3 as conda

ENV PLUGIN_SERVER=plugins.nanome.ai

COPY . /app
WORKDIR /app

RUN pip install nanome
RUN conda install -c omnia/label/cuda80 openmm
RUN conda install -c omnia pdbfixer

CMD python -m nanome_molecular_dynamics.MDSimulation -a ${PLUGIN_SERVER}

#/opt/conda/pkgs/openmm-7.4.1-py37_cuda80_rc_1/lib/python3.7/site-packages/simtk/openmm/app/data/