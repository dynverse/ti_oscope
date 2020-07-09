FROM dynverse/dynwrap_latest:v0.1.0

ARG GITHUB_PAT

RUN R -e 'devtools::install_bioc("Oscope")'

COPY definition.yml run.R example.sh /code/

ENTRYPOINT ["/code/run.R"]
