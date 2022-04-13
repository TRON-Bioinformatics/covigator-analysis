FROM jupyter/datascience-notebook:python-3.8.8


USER root


RUN python -m pip install --upgrade pip
COPY requirements.txt .
RUN pip install -r requirements.txt

RUN jupyter labextension install jupyterlab-plotly --no-build \
   && jupyter labextension install @jupyterlab/toc-extension \
   && jupyter labextension install @jupyterlab/git \
   && jupyter lab build -y \
   && jupyter lab clean -y

#RUN pip install jupyterthemes
#RUN python -m pip install --upgrade jupyterthemes
#RUN pip install jupyterlab_widgets
#RUN pip install jupyter_contrib_nbextensions
#RUN jupyter contrib nbextension install --user

# enable the Nbextensions
#RUN jupyter nbextension enable contrib_nbextensions_help_item/main
#RUN jupyter nbextension enable autosavetime/main
#RUN jupyter nbextension enable codefolding/main
#RUN jupyter nbextension enable code_font_size/code_font_size
#RUN jupyter nbextension enable code_prettify/code_prettify
#RUN jupyter nbextension enable collapsible_headings/main
#RUN jupyter nbextension enable comment-uncomment/main
#RUN jupyter nbextension enable equation-numbering/main
#RUN jupyter nbextension enable execute_time/ExecuteTime 
#RUN jupyter nbextension enable gist_it/main 
#RUN jupyter nbextension enable hide_input/main 
#RUN jupyter nbextension enable spellchecker/main
#RUN jupyter nbextension enable toc2/main
#RUN jupyter nbextension enable toggle_all_line_numbers/main


#RUN conda install --quiet --yes \
#	'conda-forge::upsetplot=0.6.0' \
#	'conda-forge::xmltodict=0.12.0' \
#	'conda-forge::biopython=1.79' \
#	'conda-forge::pybiomart=0.2.0' \
#	'conda-forge::plotly=4.14.3' \
#	'conda-forge::logzero=1.7.0'&& \
#	conda clean --all -f -y

