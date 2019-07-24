FROM rocker/shiny:latest

LABEL Description="Support for rapid design of primers for amino acid exchanges and saturation mutagenesis by Golden Gate cloning."

RUN apt-get -y update && apt-get -y install libssl-dev libxml2-dev
RUN apt-get -y --no-install-recommends install lmodern
WORKDIR /tmp
RUN wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
RUN tar -xvf install-tl-unx.tar.gz
RUN rm -rf install-tl-unx.tar.gz
ADD tex.profile .
RUN install-tl*/install-tl --profile=tex.profile
RUN wget https://github.com/jgm/pandoc/releases/download/2.7.3/pandoc-2.7.3-1-amd64.deb
ENV PATH "/usr/local/texlive/2019/bin/x86_64-linux:${PATH}"
RUN dpkg -i pandoc-2.7.3-1-amd64.deb
RUN tlmgr install dnaseq fira titlesec titling

RUN echo "PATH=/usr/local/texlive/2019/bin/x86_64-linux:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" >> /etc/environment

WORKDIR /
ADD binder/install.R /tmp
RUN R -e "source('/tmp/install.R')"

ADD GoldenMutagenesisWeb /srv/shiny-server/