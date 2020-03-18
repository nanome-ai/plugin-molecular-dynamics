FROM nvidia/cuda AS cuda
# FROM nvidia/opencl as opencl
FROM continuumio/miniconda3 as conda

ENV ARGS=''

RUN conda install -c omnia/label/cuda80 openmm
RUN conda install -c omnia pdbfixer

ARG CACHEBUST
RUN pip install nanome

COPY . /app
WORKDIR /app

CMD python run.py ${ARGS}

#/opt/conda/pkgs/openmm-7.4.1-py37_cuda80_rc_1/lib/python3.7/site-packages/simtk/openmm/app/data/
