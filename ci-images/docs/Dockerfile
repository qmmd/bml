FROM ubuntu:xenial

COPY cmake.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys DBA92F17B25AD78F9F2D9F713DEC686D130FF5E4

RUN apt-get update
RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime
RUN apt-get install --assume-yes tzdata
RUN DEBIAN_FRONTEND=noninteractive dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get install --assume-yes --no-install-recommends \
        apt-utils \
        build-essential \
        cmake cmake-data \
        doxygen \
        gfortran \
        git-core \
        libblas-dev liblapack-dev \
        openssh-client \
        texlive
