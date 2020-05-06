FROM ubuntu:xenial
RUN apt-get update
RUN apt-get install --assume-yes --no-install-recommends apt-transport-https ca-certificates wget
COPY clang.list /etc/apt/sources.list.d
COPY cmake.list /etc/apt/sources.list.d
COPY toolchain.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 60C317803A41BA51845E371A1E9377A2BA9EF27F
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys DBA92F17B25AD78F9F2D9F713DEC686D130FF5E4
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add -
RUN apt-get update
RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime
RUN apt-get install --assume-yes tzdata
RUN DEBIAN_FRONTEND=noninteractive dpkg-reconfigure --frontend noninteractive tzdata
RUN apt-get install --assume-yes --no-install-recommends apt-utils
RUN apt-get install --assume-yes --no-install-recommends build-essential
RUN apt-get install --assume-yes --no-install-recommends git-core openssh-client ca-certificates
RUN apt-get install --assume-yes --no-install-recommends cmake cmake-data libblas-dev liblapack-dev
RUN apt-get install --assume-yes --no-install-recommends gcc-4.7 g++-4.7 gfortran-4.7
RUN apt-get install --assume-yes --no-install-recommends gcc-9 g++-9 gfortran-9
RUN apt-get install --assume-yes --no-install-recommends clang-9 libomp-9-dev
RUN apt-get install --assume-yes --no-install-recommends python
RUN apt-get install --assume-yes --no-install-recommends valgrind
