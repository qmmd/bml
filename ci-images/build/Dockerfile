FROM ubuntu:xenial
COPY toolchain.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 60C317803A41BA51845E371A1E9377A2BA9EF27F
RUN apt-get update
RUN apt-get install --assume-yes --no-install-recommends apt-utils
RUN apt-get install --assume-yes --no-install-recommends build-essential
RUN apt-get install --assume-yes --no-install-recommends git-core openssh-client
RUN apt-get install --assume-yes --no-install-recommends cmake cmake-data libblas-dev liblapack-dev
RUN apt-get install --assume-yes --no-install-recommends gcc-4.7 g++-4.7 gfortran-4.7
RUN apt-get install --assume-yes --no-install-recommends gcc-9 g++-9 gfortran-9
RUN apt-get install --assume-yes --no-install-recommends python
RUN apt-get install --assume-yes --no-install-recommends valgrind
