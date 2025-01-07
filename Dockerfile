FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy the environment YAML file
COPY scmocha.prod.v4.yaml /app/scmocha.prod.v4.yaml

# Create the conda environment
RUN conda env create -f scmocha.prod.v4.yaml

# Activate the environment
SHELL ["conda", "run", "-n", "scmocha", "/bin/bash", "-c"]

# Install pip dependencies
RUN pip install \
  appdirs==1.4.4 \
  attrs==23.1.0 \
  biopython==1.81 \
  certifi==2023.7.22 \
  charset-normalizer==3.2.0 \
  click==8.1.6 \
  configargparse==1.7 \
  connection-pool==0.0.3 \
  cython==3.0.0 \
  datrie==0.8.2 \
  docutils==0.20.1 \
  dpath==2.1.6 \
  fastjsonschema==2.18.0 \
  future==0.18.3 \
  gitdb==4.0.10 \
  gitpython==3.1.32 \
  humanfriendly==10.0 \
  idna==3.4 \
  iniconfig==2.0.0 \
  jinja2==3.1.2 \
  jsonschema==4.19.0 \
  jsonschema-specifications==2023.7.1 \
  jupyter-core==5.3.1 \
  markupsafe==2.1.3 \
  mgatk==0.6.9 \
  nbformat==5.9.2 \
  numpy==1.23.5 \
  optparse-pretty==0.1.1 \
  pandas==2.0.3 \
  plac==1.3.5 \
  platformdirs==3.10.0 \
  pluggy==1.2.0 \
  psutil==5.9.5 \
  pulp==2.7.0 \
  pybktree==1.1 \
  pysam==0.21.0 \
  pytest==7.4.0 \
  pytz==2023.3 \
  referencing==0.30.2 \
  regex==2023.8.8 \
  requests==2.31.0 \
  reretry==0.11.8 \
  rpds-py==0.9.2 \
  ruamel-yaml==0.17.32 \
  ruamel-yaml-clib==0.2.7 \
  scipy==1.11.3 \
  sinto==0.10.0 \
  smart-open==6.3.0 \
  smmap==5.0.0 \
  snakemake==7.32.2 \
  stopit==1.1.2 \
  svgwrite==1.4.3 \
  tabulate==0.9.0 \
  throttler==1.2.2 \
  toposort==1.10 \
  traitlets==5.9.0 \
  tree==0.2.4 \
  umi-tools==1.1.4 \
  urllib3==2.0.4 \
  wrapt==1.15.0 \
  yte==1.5.1

# Set the entrypoint
ENTRYPOINT ["conda", "run", "-n", "scmocha", "/bin/bash", "-c"]

COPY . /opt/scMOCHA
ENV PATH /opt/scMOCHA/bin:$PATH

WORKDIR /scMOCHA