FROM rocker/shiny:latest

LABEL Description="Support for rapid design of primers for amino acid exchanges and saturation mutagenesis by Golden Gate cloning."

RUN apt-get -y update && apt-get -y install libssl-dev libxml2-dev

ADD binder/install.R /tmp
RUN R -e "source('/tmp/install.R')"

ADD GoldenMutagenesisWeb /srv/shiny-server/
ADD www /srv/shiny-server/

