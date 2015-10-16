FROM fedora:21
MAINTAINER Nuno Fonseca email: nuno.fonseca at gmail.com

# Update the image with the latest packages (recommended)
# and install missing packages
RUN sed -i "s/failovermethod=priority/failovermethod=roundrobin/" /etc/yum.repos.d/fedora.repo  && yum install -y zlib-devel python-devel bzip2-devel python python-pip R R-devel unzip make cmake git wget  tar xorg-x11-server-Xvfb && yum clean all 

# 
WORKDIR /opt

# Install and clean in a single layer
RUN git clone https://github.com/Functional-Genomics/eQTL.git  eqtl_clone  && cd /opt/eqtl_clone && ./scripts/install.sh -a /eqtl_install -s . && cd / && rm -rf /eqtl_install/tmp /opt/eqtl_clone

WORKDIR /
RUN echo source /eql_install/eqtl_setup.sh >> ~/.bash_profile && echo source /eql_install/eqtl_setup.sh >> ~/.bashrc

RUN echo '#!/usr/bin/env bash' > /usr/bin/eqtl_pipeline
RUN echo 'source /eqtl_install/eqtl_setup.sh' >> /usr/bin/eqtl_pipeline
RUN echo '/eqtl_install/scripts/eqtl_pipeline "$@"' >> /usr/bin/eqtl_pipeline
RUN chmod u+x /usr/bin/eqtl_pipeline

#ENTRYPOINT ["eql_pipeline"]

