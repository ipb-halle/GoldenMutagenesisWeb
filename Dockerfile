FROM rocker/shiny:latest

LABEL Description="Support for rapid design of primers for amino acid exchanges and saturation mutagenesis by Golden Gate cloning."

RUN apt-get -y update && apt-get -y install libssl-dev libxml2-dev pandoc
RUN apt-get -y --no-install-recommends install texlive-base texlive-latex-extra texlive-xetex texlive-fonts-extra texlive-fonts-recommended lmodern

ADD binder/install.R /tmp
RUN R -e "source('/tmp/install.R')"

ADD GoldenMutagenesisWeb /srv/shiny-server/