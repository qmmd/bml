FROM ubuntu:xenial

COPY cmake.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys DBA92F17B25AD78F9F2D9F713DEC686D130FF5E4

RUN apt-get update
RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime
RUN apt-get install --assume-yes tzdata
RUN DEBIAN_FRONTEND=noninteractive dpkg-reconfigure --frontend noninteractive tzdata
RUN apt-get install --assume-yes --no-install-recommends apt-utils
RUN apt-get install --assume-yes --no-install-recommends ca-certificates
RUN apt-get install --assume-yes --no-install-recommends build-essential
RUN apt-get install --assume-yes --no-install-recommends git-core openssh-client
RUN apt-get install --assume-yes --no-install-recommends cmake cmake-data libblas-dev liblapack-dev
RUN apt-get install --assume-yes --no-install-recommends gfortran
RUN apt-get install --assume-yes --no-install-recommends doxygen
RUN apt-get install --assume-yes --no-install-recommends texlive
