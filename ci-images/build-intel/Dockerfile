FROM ubuntu:xenial
WORKDIR /root
COPY parallel_studio_xe_2020_professional_edition_online.tgz /opt
RUN apt-get update
RUN apt-get install --assume-yes --no-install-recommends apt-utils ca-certificates
RUN apt-get install --assume-yes --no-install-recommends cpio
RUN apt-get install --assume-yes --no-install-recommends build-essential
RUN apt-get install --assume-yes --no-install-recommends linux-headers
RUN tar xvf /opt/parallel_studio_xe_2020_professional_edition_online.tgz
