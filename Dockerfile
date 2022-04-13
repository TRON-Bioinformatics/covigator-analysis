FROM jupyter/datascience-notebook:python-3.8.8


USER root


# installs python dependencies
RUN python -m pip install --upgrade pip
COPY requirements.txt .
RUN pip install -r requirements.txt

# installs some extensions to Jupyter
RUN jupyter labextension install jupyterlab-plotly --no-build \
   && jupyter labextension install @jupyterlab/toc-extension \
   && jupyter labextension install @jupyterlab/git \
   && jupyter lab build -y \
   && jupyter lab clean -y

