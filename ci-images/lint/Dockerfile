FROM ubuntu:xenial

COPY emacs.list /etc/apt/sources.list.d
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 873503A090750CDAEB0754D93FF0E01EEAAFC9CD

RUN apt-get update
RUN ln -fs /usr/share/zoneinfo/UTC /etc/localtime
RUN apt-get install --assume-yes tzdata
RUN DEBIAN_FRONTEND=noninteractive dpkg-reconfigure --frontend noninteractive tzdata
RUN apt-get install --assume-yes --no-install-recommends apt-utils
RUN apt-get install --assume-yes --no-install-recommends ca-certificates
RUN apt-get install --assume-yes --no-install-recommends git-core openssh-client
RUN apt-get install --assume-yes --no-install-recommends bundler
RUN apt-get install --assume-yes --no-install-recommends emacs26
RUN apt-get install --assume-yes --no-install-recommends indent
RUN apt-get install --assume-yes --no-install-recommends python-pip
RUN apt-get install --assume-yes --no-install-recommends python-setuptools
RUN apt-get install --assume-yes --no-install-recommends python-wheel
RUN pip install proselint
