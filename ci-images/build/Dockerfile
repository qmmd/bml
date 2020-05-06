FROM ubuntu:xenial

RUN apt-get update
RUN apt-get install --assume-yes --no-install-recommends apt-transport-https ca-certificates wget

COPY clang.list /etc/apt/sources.list.d
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add -
COPY cmake.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys DBA92F17B25AD78F9F2D9F713DEC686D130FF5E4
COPY toolchain.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 60C317803A41BA51845E371A1E9377A2BA9EF27F

RUN apt-get update
RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime
RUN apt-get install --assume-yes tzdata
RUN DEBIAN_FRONTEND=noninteractive dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get install --assume-yes --no-install-recommends \
        apt-utils \
        build-essential \
        clang-9 llvm-9-dev libomp-9-dev \
        cmake cmake-data \
        gcc-4.7 g++-4.7 gfortran-4.7 \
        gcc-9 g++-9 gfortran-9 \
        git-core openssh-client \
        libblas-dev liblapack-dev \
        python \
        valgrind
